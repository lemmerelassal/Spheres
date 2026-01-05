// PDBViewer.cpp — UE 5.6 compatible version

#include "PDBViewer.h"
#include "PDBCameraComponent.h"
#include "HttpModule.h"
#include "Interfaces/IHttpResponse.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Engine/StaticMesh.h"
#include "Engine/World.h"
#include "Components/StaticMeshComponent.h"
#include "Components/SceneComponent.h"
#include "UObject/ConstructorHelpers.h"
#include "Misc/FileHelper.h"
#include "HAL/PlatformFilemanager.h"
#include "Misc/Paths.h"
#include "IDesktopPlatform.h"
#include "DesktopPlatformModule.h"
#include "Interfaces/IMainFrameModule.h"

namespace PDB
{
    constexpr float SCALE = 50.0f;
    constexpr float SPHERE_SIZE = 0.5f;
    constexpr float CYLINDER_SIZE = 0.1f;
    constexpr float BOND_OFFSET = 8.0f;
}

APDBViewer::APDBViewer()
{
    PrimaryActorTick.bCanEverTick = false;
    RootComponent = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));

    static ConstructorHelpers::FObjectFinder<UStaticMesh> Sphere(TEXT("/Engine/BasicShapes/Sphere"));
    static ConstructorHelpers::FObjectFinder<UStaticMesh> Cylinder(TEXT("/Engine/BasicShapes/Cylinder"));
    static ConstructorHelpers::FObjectFinder<UMaterial> Mat(TEXT("/Engine/BasicShapes/BasicShapeMaterial"));

    SphereMeshAsset = Sphere.Object;
    CylinderMeshAsset = Cylinder.Object;
    SphereMaterialAsset = Mat.Object;
}

void APDBViewer::BeginPlay()
{
    Super::BeginPlay();
    FetchAndDisplayStructure(TEXT("5ENB"));

    if (auto* Cam = GetWorld()->SpawnActor<APDBCameraComponent>(APDBCameraComponent::StaticClass(), GetActorLocation(), FRotator::ZeroRotator))
        Cam->SetTargetActor(this);
}

void APDBViewer::FetchAndDisplayStructure(const FString& ID)
{
    CurrentStructureID = ID;
    FString URL = FString::Printf(TEXT("https://files.rcsb.org/download/%s.pdb"), *ID);
    FetchFileAsync(URL, [this, ID](bool bOK, const FString& Content) {
        if (bOK) ParsePDB(Content);
        else FetchFileAsync(FString::Printf(TEXT("https://files.rcsb.org/download/%s.cif"), *ID),
            [this](bool bOK2, const FString& C) { if (bOK2) ParseMMCIF(C); });
    });
}

void APDBViewer::FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> CB)
{
    auto Req = FHttpModule::Get().CreateRequest();
    Req->SetURL(URL);
    Req->SetVerb(TEXT("GET"));
    Req->OnProcessRequestComplete().BindLambda([CB](FHttpRequestPtr R, FHttpResponsePtr Resp, bool bOK) {
        CB(bOK && Resp.IsValid() && Resp->GetResponseCode() == 200, bOK ? Resp->GetContentAsString() : TEXT(""));
    });
    Req->ProcessRequest();
}

void APDBViewer::ParsePDB(const FString& Content)
{
    CurrentPDBContent = Content;
    ClearResidueMap();

    TArray<FString> Lines;
    Content.ParseIntoArrayLines(Lines);

    TMap<FString, TMap<FString, FVector>> ResAtoms;
    TMap<FString, FResidueMetadata> ResMeta;
    FString FirstChain;

    for (const auto& L : Lines)
    {
        if (L.Len() < 80 || !L.StartsWith(TEXT("HETATM"))) continue;

        FString Chain = L.Mid(21, 1);
        if (FirstChain.IsEmpty()) FirstChain = Chain;
        else if (Chain != FirstChain) continue;

        FString Key = FString::Printf(TEXT("%s_%s_%s"),
            *L.Mid(17, 3).TrimStartAndEnd(),
            *L.Mid(22, 4).TrimStartAndEnd(),
            *Chain);

        ResAtoms.FindOrAdd(Key).Add(L.Mid(12, 4).TrimStartAndEnd(),
            FVector(FCString::Atof(*L.Mid(30, 8)),
                FCString::Atof(*L.Mid(38, 8)),
                FCString::Atof(*L.Mid(46, 8))) * PDB::SCALE);

        if (!ResMeta.Contains(Key))
        {
            auto& M = ResMeta.Add(Key);
            M.ResidueName = L.Mid(17, 3).TrimStartAndEnd();
            M.ResidueSeq = L.Mid(22, 4).TrimStartAndEnd();
            M.Chain = Chain;
            M.RecordType = TEXT("HETATM");
        }
    }

    CreateResiduesFromAtomData(ResAtoms, ResMeta);
    OnResiduesLoaded.Broadcast();
}

void APDBViewer::ParseMMCIF(const FString& Content)
{
    CurrentPDBContent = Content;
    ClearResidueMap();

    TArray<FString> Lines;
    Content.ParseIntoArrayLines(Lines);

    TArray<FString> Hdrs;
    TArray<TArray<FString>> AtomTab;
    int32 XI = -1, YI = -1, ZI = -1, RI = -1, AI = -1, GI = -1, CI = -1, SI = -1;
    bool bLoop = false;

    for (const auto& L : Lines)
    {
        if (L.StartsWith(TEXT("loop_"))) { bLoop = true; Hdrs.Empty(); continue; }
        if (bLoop && L.StartsWith(TEXT("_atom_site.")))
        {
            int32 I = Hdrs.Add(L);
            if (L.Contains(TEXT("Cartn_x"))) XI = I;
            else if (L.Contains(TEXT("Cartn_y"))) YI = I;
            else if (L.Contains(TEXT("Cartn_z"))) ZI = I;
            else if (L.Contains(TEXT("label_comp_id"))) RI = I;
            else if (L.Contains(TEXT("label_atom_id"))) AI = I;
            else if (L.Contains(TEXT("group_PDB"))) GI = I;
            else if (L.Contains(TEXT("label_asym_id"))) CI = I;
            else if (L.Contains(TEXT("label_seq_id")) || L.Contains(TEXT("auth_seq_id"))) SI = I;
            continue;
        }
        if (bLoop && !L.StartsWith(TEXT("_")))
        {
            TArray<FString> T;
            L.ParseIntoArrayWS(T);
            if (T.Num() > FMath::Max3(XI, YI, ZI)) AtomTab.Add(MoveTemp(T));
        }
    }

    TMap<FString, TMap<FString, FVector>> ResAtoms;
    TMap<FString, FResidueMetadata> ResMeta;
    FString FirstChain;

    for (const auto& R : AtomTab)
    {
        if (XI < 0 || YI < 0 || ZI < 0 || RI < 0 || AI < 0) continue;
        if (GI >= 0 && R.IsValidIndex(GI) && R[GI] != TEXT("ATOM")) continue;

        FString Chain = (CI >= 0 && R.IsValidIndex(CI)) ? R[CI] : TEXT("");
        if (FirstChain.IsEmpty()) FirstChain = Chain;
        else if (Chain != FirstChain) continue;

        FString Seq = (SI >= 0 && R.IsValidIndex(SI)) ? R[SI] : TEXT("0");
        FString Key = FString::Printf(TEXT("%s_%s_%s"), *R[RI], *Seq, *Chain);

        ResAtoms.FindOrAdd(Key).Add(R[AI],
            FVector(FCString::Atof(*R[XI]),
                FCString::Atof(*R[YI]),
                FCString::Atof(*R[ZI])) * PDB::SCALE);

        if (!ResMeta.Contains(Key))
        {
            auto& M = ResMeta.Add(Key);
            M.ResidueName = R[RI];
            M.ResidueSeq = Seq;
            M.Chain = Chain;
            M.RecordType = TEXT("ATOM");
        }
    }

    CreateResiduesFromAtomData(ResAtoms, ResMeta);
    OnResiduesLoaded.Broadcast();
}

void APDBViewer::CreateResiduesFromAtomData(const TMap<FString, TMap<FString, FVector>>& ResAtoms, const TMap<FString, FResidueMetadata>& Meta)
{
    for (const auto& P : ResAtoms)
    {
        const auto* M = Meta.Find(P.Key);
        if (!M) continue;

        auto* Info = new FResidueInfo();
        Info->ResidueName = M->ResidueName;
        Info->ResidueSeq = M->ResidueSeq;
        Info->Chain = M->Chain;
        Info->RecordType = M->RecordType;
        Info->bIsVisible = true;

        for (const auto& A : P.Value)
            DrawSphere(A.Value.X, A.Value.Y, A.Value.Z, GetElementColor(A.Key.Left(1)), GetRootComponent(), Info->AtomMeshes);

        ResidueMap.Add(P.Key, Info);
        FetchLigandBondsForResidue(P.Key, M->ResidueName, P.Value);
    }
}

void APDBViewer::FetchLigandBondsForResidue(const FString& Key, const FString& Name, const TMap<FString, FVector>& Pos)
{
    FetchFileAsync(FString::Printf(TEXT("https://files.rcsb.org/ligands/download/%s.cif"), *Name.ToUpper()),
        [this, Key, Pos](bool bOK, const FString& C) {
            if (bOK && ResidueMap.Contains(Key))
                if (auto** Info = ResidueMap.Find(Key))
                    if (*Info) ParseLigandCIFForResidue(C, Pos, *Info);
        });
}

void APDBViewer::ParseLigandCIFForResidue(const FString& Content, const TMap<FString, FVector>& Pos, FResidueInfo* Info)
{
    if (!Info) return;

    TMap<FString, FVector> NormPos;
    for (const auto& P : Pos)
    {
        FString K;
        for (const TCHAR C : P.Key) if (FChar::IsAlnum(C)) K.AppendChar(FChar::ToUpper(C));
        if (!K.IsEmpty()) NormPos.Add(K, P.Value);
    }

    TArray<FString> Lines;
    Content.ParseIntoArrayLines(Lines);

    TArray<FString> Hdrs;
    int32 A1 = -1, A2 = -1, BO = -1;
    bool bLoop = false;

    for (const auto& L : Lines)
    {
        if (L.StartsWith(TEXT("loop_"))) { bLoop = true; Hdrs.Empty(); A1 = A2 = BO = -1; continue; }
        if (bLoop && L.StartsWith(TEXT("_")))
        {
            int32 I = Hdrs.Add(L);
            FString Lo = L.ToLower();
            if (Lo.Contains(TEXT("atom_id_1")) || Lo.Contains(TEXT("atom_1"))) A1 = I;
            else if (Lo.Contains(TEXT("atom_id_2")) || Lo.Contains(TEXT("atom_2"))) A2 = I;
            else if (Lo.Contains(TEXT("value_order")) || Lo.Contains(TEXT("bond_order"))) BO = I;
            continue;
        }
        if (bLoop && !L.StartsWith(TEXT("_")) && !L.StartsWith(TEXT("data_")))
        {
            if (A1 < 0 || A2 < 0) continue;
            TArray<FString> T;
            L.ParseIntoArrayWS(T);
            if (T.Num() <= FMath::Max(A1, A2)) continue;

            FString ID1 = NormalizeAtomID(T[A1]);
            FString ID2 = NormalizeAtomID(T[A2]);
            int32 Ord = ParseBondOrder(BO >= 0 && T.IsValidIndex(BO) ? T[BO] : TEXT("1"));

            const auto* P1 = NormPos.Find(ID1);
            const auto* P2 = NormPos.Find(ID2);
            if (P1 && P2) DrawBond(*P1, *P2, Ord, FLinearColor::Gray, GetRootComponent(), Info->BondMeshes);
        }
        else if (bLoop && L.StartsWith(TEXT("data_"))) break;
    }
}

FString APDBViewer::NormalizeAtomID(const FString& In) const
{
    FString Out;
    for (const TCHAR C : In) if (FChar::IsAlnum(C)) Out.AppendChar(FChar::ToUpper(C));
    return Out;
}

int32 APDBViewer::ParseBondOrder(const FString& S) const
{
    if (S.Len() == 1 && S[0] >= '1' && S[0] <= '3') return S[0] - '0';
    FString L = S.ToLower();
    if (L.Contains(TEXT("ar"))) return 1;
    if (L.Contains(TEXT("doub")) || L.Equals(TEXT("d"))) return 2;
    if (L.Contains(TEXT("trip")) || L.Equals(TEXT("t"))) return 3;
    return FMath::Clamp(FCString::Atoi(*S), 1, 3);
}

void APDBViewer::DrawSphere(float X, float Y, float Z, const FLinearColor& Col, USceneComponent* Par, TArray<UStaticMeshComponent*>& Out)
{
    if (!SphereMeshAsset || !SphereMaterialAsset || !Par) return;

    auto* Sph = NewObject<UStaticMeshComponent>(this);
    Sph->SetStaticMesh(SphereMeshAsset);
    Sph->SetWorldScale3D(FVector(PDB::SPHERE_SIZE));
    Sph->SetCollisionEnabled(ECollisionEnabled::NoCollision);

    auto* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, Sph);
    Mat->SetVectorParameterValue(TEXT("Color"), Col);
    Mat->SetScalarParameterValue(TEXT("EmissiveIntensity"), 5.0f);
    Sph->SetMaterial(0, Mat);

    Sph->AttachToComponent(Par, FAttachmentTransformRules::KeepWorldTransform);
    Sph->SetRelativeLocation(Par->GetComponentTransform().InverseTransformPosition(FVector(X, Y, Z)));
    Sph->RegisterComponent();

    Out.Add(Sph);
    AllAtomMeshes.Add(Sph);
}

void APDBViewer::DrawBond(const FVector& S, const FVector& E, int32 Ord, const FLinearColor& Col, USceneComponent* Par, TArray<UStaticMeshComponent*>& Out)
{
    if (!CylinderMeshAsset || !SphereMaterialAsset || !Par) return;

    FVector V = E - S;
    float Len = V.Size();
    if (Len < KINDA_SMALL_NUMBER) return;

    FVector Mid = S + V * 0.5f;
    FRotator Rot = FRotationMatrix::MakeFromZ(V).Rotator();
    float Scale = Len / 100.0f;

    auto MakeCyl = [&](const FVector& Pos)
    {
        auto* Cyl = NewObject<UStaticMeshComponent>(this);
        Cyl->SetStaticMesh(CylinderMeshAsset);
        Cyl->SetWorldLocation(Pos);
        Cyl->SetWorldRotation(Rot);
        Cyl->SetWorldScale3D(FVector(PDB::CYLINDER_SIZE, PDB::CYLINDER_SIZE, Scale));
        Cyl->SetCollisionEnabled(ECollisionEnabled::NoCollision);
        auto* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, Cyl);
        Mat->SetVectorParameterValue(TEXT("Color"), Col);
        Cyl->SetMaterial(0, Mat);
        Cyl->AttachToComponent(Par, FAttachmentTransformRules::KeepWorldTransform);
        Cyl->RegisterComponent();
        Out.Add(Cyl);
        AllBondMeshes.Add(Cyl);
    };

    if (Ord <= 1) { MakeCyl(Mid); return; }

    FVector Perp = FVector::CrossProduct(V.GetSafeNormal(), FVector::UpVector);
    if (Perp.SizeSquared() < KINDA_SMALL_NUMBER) Perp = FVector::CrossProduct(V.GetSafeNormal(), FVector::RightVector);
    Perp.Normalize();
    FVector Off = Perp * PDB::BOND_OFFSET;

    if (Ord == 2) { MakeCyl(Mid + Off); MakeCyl(Mid - Off); }
    else if (Ord == 3) { MakeCyl(Mid); MakeCyl(Mid + Off); MakeCyl(Mid - Off); }
    else MakeCyl(Mid);
}

FLinearColor APDBViewer::GetElementColor(const FString& E) const
{
    static const TMap<FString, FLinearColor> Colors = {
        {TEXT("C"), FLinearColor(0.1f,0.1f,0.1f)}, {TEXT("O"), FLinearColor::Red}, {TEXT("H"), FLinearColor::White},
        {TEXT("D"), FLinearColor::White}, {TEXT("N"), FLinearColor::Blue}, {TEXT("S"), FLinearColor::Yellow},
        {TEXT("CL"), FLinearColor(0,1,0)}, {TEXT("P"), FLinearColor(1,0.5f,0)}, {TEXT("F"), FLinearColor(0,1,0)},
        {TEXT("BR"), FLinearColor(0.6f,0.2f,0.2f)}, {TEXT("I"), FLinearColor(0.4f,0,0.8f)},
        {TEXT("FE"), FLinearColor(0.8f,0.4f,0)}, {TEXT("MG"), FLinearColor(0,0.8f,0)},
        {TEXT("ZN"), FLinearColor(0.5f,0.5f,0.5f)}, {TEXT("CA"), FLinearColor(0.2f,0.6f,1)},
        {TEXT("NA"), FLinearColor(0,0,1)}, {TEXT("K"), FLinearColor(0.5f,0,1)}, {TEXT("CU"), FLinearColor(1,0.5f,0)}, {TEXT("B"), FLinearColor(1,0.7f,0.7f)}
    };
    const auto* C = Colors.Find(E.ToUpper());
    return C ? *C : FLinearColor::Gray;
}

void APDBViewer::ClearResidueMap()
{
    for (auto& P : ResidueMap) delete P.Value;
    ResidueMap.Empty();
}

void APDBViewer::ClearCurrentStructure()
{
    ClearResidueMap();
    for (auto* M : AllAtomMeshes) if (M && IsValid(M)) M->DestroyComponent();
    for (auto* M : AllBondMeshes) if (M && IsValid(M)) M->DestroyComponent();
    AllAtomMeshes.Empty();
    AllBondMeshes.Empty();
}

void APDBViewer::ToggleResidueVisibility(const FString& Key)
{
    auto** InfoPtr = ResidueMap.Find(Key);
    if (!InfoPtr || !*InfoPtr) return;
    auto* Info = *InfoPtr;
    Info->bIsVisible = !Info->bIsVisible;
    for (auto* M : Info->AtomMeshes) if (M) M->SetVisibility(Info->bIsVisible);
    for (auto* M : Info->BondMeshes) if (M) M->SetVisibility(Info->bIsVisible);
}

TArray<FString> APDBViewer::GetResidueList() const
{
    TArray<FString> List;
    ResidueMap.GetKeys(List);
    List.Sort([](const FString& A, const FString& B) {
        int32 U1, U2;
        if (!A.FindChar('_', U1) || !B.FindChar('_', U2))
            return A < B;

        int32 U3 = A.Find(TEXT("_"), ESearchCase::IgnoreCase, ESearchDir::FromStart, U1 + 1);
        int32 U4 = B.Find(TEXT("_"), ESearchCase::IgnoreCase, ESearchDir::FromStart, U2 + 1);

        if (U3 != INDEX_NONE && U4 != INDEX_NONE)
        {
            int32 NumA = FCString::Atoi(*A.Mid(U1 + 1, U3 - U1 - 1));
            int32 NumB = FCString::Atoi(*B.Mid(U2 + 1, U4 - U2 - 1));
            return NumA < NumB;
        }
        return A < B;
    });
    return List;
}

TArray<FString> APDBViewer::GetLigandList() const { return GetResidueList(); }

FString APDBViewer::GetResidueDisplayName(const FString& Key) const
{
    const auto* P = ResidueMap.Find(Key);
    return (P && *P)
        ? FString::Printf(TEXT("%s %s"), *(*P)->ResidueName, *(*P)->ResidueSeq)
        : Key;
}

void APDBViewer::SaveStructureToFile(const FString& Path)
{
    if (CurrentPDBContent.IsEmpty())
    {
        UE_LOG(LogTemp, Warning, TEXT("No data to save"));
        return;
    }

    if (FFileHelper::SaveStringToFile(CurrentPDBContent, *Path))
    {
        UE_LOG(LogTemp, Log, TEXT("Saved: %s"), *Path);
    }
    else
    {
        UE_LOG(LogTemp, Error, TEXT("Failed: %s"), *Path);
    }
}

void APDBViewer::LoadStructureFromFile(const FString& Path)
{
    FString Content;
    if (!FFileHelper::LoadFileToString(Content, *Path))
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to load: %s"), *Path);
        return;
    }

    ClearCurrentStructure();
    CurrentStructureID = FPaths::GetBaseFilename(Path);
    CurrentPDBContent = MoveTemp(Content);

    FString Ext = FPaths::GetExtension(Path).ToLower();
    if (Ext == TEXT("pdb")) ParsePDB(CurrentPDBContent);
    else if (Ext == TEXT("cif")) ParseMMCIF(CurrentPDBContent);
    else UE_LOG(LogTemp, Error, TEXT("Unsupported format: %s"), *Ext);
}

bool APDBViewer::ShowFileDialog(bool bSave, FString& Out)
{
    auto* DP = FDesktopPlatformModule::Get();
    if (!DP) return false;

    void* Handle = nullptr;
    if (FModuleManager::Get().IsModuleLoaded("MainFrame"))
        if (auto P = FModuleManager::LoadModuleChecked<IMainFrameModule>("MainFrame").GetParentWindow())
            if (P.IsValid() && P->GetNativeWindow().IsValid())
                Handle = P->GetNativeWindow()->GetOSWindowHandle();

    TArray<FString> Files;
    bool bOK = bSave ?
        DP->SaveFileDialog(Handle, TEXT("Save PDB"), FPaths::ProjectSavedDir(), TEXT("structure.pdb"),
            TEXT("PDB Files (*.pdb)|*.pdb|mmCIF Files (*.cif)|*.cif"), EFileDialogFlags::None, Files) :
        DP->OpenFileDialog(Handle, TEXT("Load PDB"), FPaths::ProjectSavedDir(), TEXT(""),
            TEXT("PDB Files (*.pdb)|*.pdb|mmCIF Files (*.cif)|*.cif"), EFileDialogFlags::None, Files);

    if (bOK && Files.Num() > 0) { Out = MoveTemp(Files[0]); return true; }
    return false;
}

void APDBViewer::OpenSaveDialog() { FString P; if (ShowFileDialog(true, P)) SaveStructureToFile(P); }
void APDBViewer::OpenLoadDialog() { FString P; if (ShowFileDialog(false, P)) LoadStructureFromFile(P); }
