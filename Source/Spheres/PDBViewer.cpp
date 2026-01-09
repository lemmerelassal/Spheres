// PDBViewer.cpp – UE 5.6 compatible version with TreeView and SDF support

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
#include "Components/TreeView.h"
#include "Components/ListView.h"

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
    ChainIDs.Empty();

    TArray<FString> Lines;
    Content.ParseIntoArrayLines(Lines);

    TMap<FString, TMap<FString, FVector>> ResAtoms;
    TMap<FString, FResidueMetadata> ResMeta;

    for (const auto& L : Lines)
    {
        FString RecordType = L.Mid(0, 6);
        if (L.Len() < 80  || ! (RecordType.StartsWith("ATOM") || RecordType.StartsWith("HETATM")) ) continue;

        FString Chain = L.Mid(21, 1).TrimStartAndEnd();
        if (Chain.IsEmpty()) Chain = TEXT("_"); // Handle empty chain IDs
        
        ChainIDs.Add(Chain);

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
            M.RecordType = RecordType;
        }
    }

    CreateResiduesFromAtomData(ResAtoms, ResMeta);
    OnResiduesLoaded.Broadcast();
}

void APDBViewer::ParseMMCIF(const FString& Content)
{
    CurrentPDBContent = Content;
    ClearResidueMap();
    ChainIDs.Empty();

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

    for (const auto& R : AtomTab)
    {
        if (XI < 0 || YI < 0 || ZI < 0 || RI < 0 || AI < 0) continue;
        if (GI >= 0 && R.IsValidIndex(GI) && R[GI] != TEXT("ATOM")) continue;

        FString Chain = (CI >= 0 && R.IsValidIndex(CI)) ? R[CI] : TEXT("_");
        if (Chain.IsEmpty()) Chain = TEXT("_");
        
        ChainIDs.Add(Chain);

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

void APDBViewer::ParseSDF(const FString& Content)
{
    CurrentPDBContent = Content;
    ClearLigandMap();

    TArray<FString> Lines;
    Content.ParseIntoArrayLines(Lines);

    int32 MoleculeIndex = 0;
    int32 LineIndex = 0;
    
    while (LineIndex < Lines.Num())
    {
        // Check if we have enough lines for a molecule header
        if (LineIndex + 3 >= Lines.Num())
            break;

        FString MoleculeName = Lines[LineIndex].TrimStartAndEnd();
        if (MoleculeName.IsEmpty())
            MoleculeName = FString::Printf(TEXT("MOL%d"), MoleculeIndex + 1);
        
        // Find the counts line - it's the line with the V2000 or V3000 marker
        // Usually at LineIndex + 3, but let's search for it
        int32 CountsLineIndex = -1;
        for (int32 i = LineIndex + 1; i < FMath::Min(LineIndex + 5, Lines.Num()); ++i)
        {
            if (Lines[i].Contains(TEXT("V2000")) || Lines[i].Contains(TEXT("V3000")))
            {
                CountsLineIndex = i;
                break;
            }
        }
        
        if (CountsLineIndex == -1)
        {
            UE_LOG(LogTemp, Warning, TEXT("Could not find counts line for molecule '%s'"), *MoleculeName);
            LineIndex++;
            continue;
        }
        
        FString CountsLine = Lines[CountsLineIndex];
        
        UE_LOG(LogTemp, Log, TEXT("SDF molecule '%s' at line %d"), *MoleculeName, LineIndex);
        UE_LOG(LogTemp, Log, TEXT("  Counts line (index %d): '%s' (len=%d)"), CountsLineIndex, *CountsLine, CountsLine.Len());
        
        if (CountsLine.Len() < 6)
        {
            UE_LOG(LogTemp, Warning, TEXT("Invalid counts line, skipping molecule"));
            LineIndex++;
            continue;
        }

        int32 NumAtoms = FCString::Atoi(*CountsLine.Mid(0, 3).TrimStartAndEnd());
        int32 NumBonds = FCString::Atoi(*CountsLine.Mid(3, 3).TrimStartAndEnd());

        UE_LOG(LogTemp, Log, TEXT("  Expecting %d atoms, %d bonds"), NumAtoms, NumBonds);

        // Parse atom block (starts right after counts line)
        TArray<FVector> AtomPositions;
        TArray<FString> AtomElements;
        int32 AtomStartLine = CountsLineIndex + 1;
        
        UE_LOG(LogTemp, Log, TEXT("  AtomStartLine=%d, Lines.Num()=%d"), AtomStartLine, Lines.Num());
        
        for (int32 i = 0; i < NumAtoms; ++i)
        {
            int32 CurrentLine = AtomStartLine + i;
            if (CurrentLine >= Lines.Num())
            {
                UE_LOG(LogTemp, Warning, TEXT("  Ran out of lines at atom %d (line %d >= %d)"), i, CurrentLine, Lines.Num());
                break;
            }
            
            FString Line = Lines[CurrentLine];
            if (Line.Len() < 34) 
            {
                UE_LOG(LogTemp, Warning, TEXT("  Atom line %d too short: '%s' (len=%d)"), i, *Line, Line.Len());
                continue;
            }

            // Parse coordinates
            float X = FCString::Atof(*Line.Mid(0, 10).TrimStartAndEnd());
            float Y = FCString::Atof(*Line.Mid(10, 10).TrimStartAndEnd());
            float Z = FCString::Atof(*Line.Mid(20, 10).TrimStartAndEnd());
            
            // Parse element symbol (starts at position 31)
            FString Element = Line.Mid(31, 3).TrimStartAndEnd();
            FVector Pos(X, Y, Z);
            Pos *= PDB::SCALE;
            
            AtomPositions.Add(Pos);
            AtomElements.Add(Element);
            
            if (i == 0)
            {
                UE_LOG(LogTemp, Log, TEXT("  First atom: Element='%s' Pos=(%f,%f,%f) Scaled=(%f,%f,%f)"), 
                       *Element, X, Y, Z, Pos.X, Pos.Y, Pos.Z);
            }
        }
        
        UE_LOG(LogTemp, Log, TEXT("  Parsed %d atoms (expected %d)"), AtomPositions.Num(), NumAtoms);

        // Parse bond block (starts after atom block)
        int32 BondStartLine = AtomStartLine + NumAtoms;
        TArray<TPair<int32, int32>> BondPairs;
        TArray<int32> BondOrders;

        for (int32 i = 0; i < NumBonds && (BondStartLine + i) < Lines.Num(); ++i)
        {
            FString Line = Lines[BondStartLine + i];
            if (Line.Len() < 9) continue;

            // Parse bond line
            int32 Atom1 = FCString::Atoi(*Line.Mid(0, 3).TrimStartAndEnd()) - 1; // Convert to 0-indexed
            int32 Atom2 = FCString::Atoi(*Line.Mid(3, 3).TrimStartAndEnd()) - 1;
            int32 BondType = FCString::Atoi(*Line.Mid(6, 3).TrimStartAndEnd());

            if (Atom1 >= 0 && Atom1 < AtomPositions.Num() && 
                Atom2 >= 0 && Atom2 < AtomPositions.Num())
            {
                BondPairs.Add(TPair<int32, int32>(Atom1, Atom2));
                BondOrders.Add(BondType);
            }
        }

        // Create residue for this molecule
        FString Key = FString::Printf(TEXT("%s"), *MoleculeName);
        
        auto* Info = new FLigandInfo();
        Info->LigandName = MoleculeName;
        Info->bIsVisible = true;

        // Draw all atoms
        for (int32 i = 0; i < AtomPositions.Num(); ++i)
        {
            FLinearColor Color = GetElementColor(AtomElements[i]);
            DrawSphere(AtomPositions[i].X, AtomPositions[i].Y, AtomPositions[i].Z, 
                       Color, GetRootComponent(), Info->AtomMeshes);
        }

        // Draw all bonds
        for (int32 i = 0; i < BondPairs.Num(); ++i)
        {
            int32 Atom1 = BondPairs[i].Key;
            int32 Atom2 = BondPairs[i].Value;
            int32 Order = BondOrders[i];

            // Treat aromatic bonds as single bonds for visualization
            if (Order == 4) Order = 1;
            
            DrawBond(AtomPositions[Atom1], AtomPositions[Atom2], Order, 
                     FLinearColor::Gray, GetRootComponent(), Info->BondMeshes);
        }

        UE_LOG(LogTemp, Log, TEXT("  Created %d atom meshes, %d bond meshes for '%s'"), 
               Info->AtomMeshes.Num(), Info->BondMeshes.Num(), *Key);

        LigandMap.Add(Key, Info);
        
        // Find the next molecule separator ($$$$) or end of file
        LineIndex = BondStartLine + NumBonds;
        while (LineIndex < Lines.Num())
        {
            if (Lines[LineIndex].StartsWith(TEXT("$$")))
            {
                LineIndex++; // Skip the $$$$ line
                break;
            }
            LineIndex++;
        }
        
        MoleculeIndex++;
    }
    
    OnLigandsLoaded.Broadcast();
    
    UE_LOG(LogTemp, Log, TEXT("SDF loaded successfully: %d molecules"), MoleculeIndex);
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
    ChainIDs.Empty();
}

void APDBViewer::ClearLigandMap() {
    for (auto& P : LigandMap) delete P.Value;
    LigandMap.Empty();
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

// New TreeView Functions

TArray<UPDBTreeNode*> APDBViewer::GetChainNodes()
{
    TArray<UPDBTreeNode*> Nodes;
    TArray<FString> SortedChains = ChainIDs.Array();
    SortedChains.Sort();
    
    for (const FString& ChainID : SortedChains)
    {
        FString DisplayName = ChainID == TEXT("_") 
            ? TEXT("Chain (No ID)") 
            : FString::Printf(TEXT("Chain %s"), *ChainID);
        
        UPDBTreeNode* Node = NewObject<UPDBTreeNode>(this);
        Node->Initialize(DisplayName, ChainID, true, ChainID);
        Nodes.Add(Node);
    }
    
    return Nodes;
}

TArray<UPDBTreeNode*> APDBViewer::GetResidueNodesForChain(const FString& ChainID)
{
    TArray<UPDBTreeNode*> Nodes;
    TArray<FString> Keys;
    ResidueMap.GetKeys(Keys);
    
    // Sort by residue sequence number
    Keys.Sort([](const FString& A, const FString& B) {
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
    
    for (const FString& Key : Keys)
    {
        const auto* InfoPtr = ResidueMap.Find(Key);
        if (!InfoPtr || !*InfoPtr) continue;
        
        const FResidueInfo* Info = *InfoPtr;
        if (Info->Chain != ChainID) continue;
        
        FString DisplayName = FString::Printf(TEXT("%s %s"), *Info->ResidueName, *Info->ResidueSeq);
        
        UPDBTreeNode* Node = NewObject<UPDBTreeNode>(this);
        Node->Initialize(DisplayName, Key, false, ChainID);
        Node->bIsVisible = Info->bIsVisible;
        
        Nodes.Add(Node);
    }
    
    return Nodes;
}

void APDBViewer::ToggleChainVisibility(const FString& ChainID)
{
    bool bNewVisibility = false;
    bool bFirstResidue = true;
    
    // Determine new visibility state based on first residue in chain
    for (const auto& Pair : ResidueMap)
    {
        if (Pair.Value && Pair.Value->Chain == ChainID)
        {
            if (bFirstResidue)
            {
                bNewVisibility = !Pair.Value->bIsVisible;
                bFirstResidue = false;
            }
            
            Pair.Value->bIsVisible = bNewVisibility;
            for (auto* M : Pair.Value->AtomMeshes) if (M) M->SetVisibility(bNewVisibility);
            for (auto* M : Pair.Value->BondMeshes) if (M) M->SetVisibility(bNewVisibility);
        }
    }
}

void APDBViewer::ToggleNodeVisibility(UPDBTreeNode* Node)
{
    if (!Node) return;
    
    if (Node->bIsChain)
    {
        ToggleChainVisibility(Node->ChainID);
    }
    else
    {
        ToggleResidueVisibility(Node->NodeKey);
    }
}

void APDBViewer::PopulateTreeView(UTreeView* TreeView)
{
    if (!TreeView) return;
    
    // Clear existing items
    TreeView->ClearListItems();
    
    // Get all chain nodes
    TArray<UPDBTreeNode*> ChainNodes = GetChainNodes();
    
    // Add each chain to the tree view
    for (UPDBTreeNode* ChainNode : ChainNodes)
    {
        TreeView->AddItem(ChainNode);
    }
    
    // Request refresh
    TreeView->RequestRefresh();
}

TArray<UObject*> APDBViewer::GetChildrenForNode(UPDBTreeNode* Node)
{
    TArray<UObject*> ChildNodes;
    
    if (!Node)
        return ChildNodes;
    
    if (Node->bIsChain)
    {
        // Get residues for this chain
        TArray<UPDBTreeNode*> Residues = GetResidueNodesForChain(Node->ChainID);
        
        // Convert to UObject array
        for (UPDBTreeNode* Residue : Residues)
        {
            ChildNodes.Add(Residue);
        }
    }
    // else: residues have no children, return empty array
    
    return ChildNodes;
}

// Legacy Functions (kept for backwards compatibility)

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


FString APDBViewer::GetLigandDisplayName(const FString& Key) const
{
    const auto* P = LigandMap.Find(Key);
    return (P && *P)
        ? FString::Printf(TEXT("%s"), *(*P)->LigandName)
        : Key;
}

// File I/O Functions

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
    if (Ext == TEXT("pdb")) 
        ParsePDB(CurrentPDBContent);
    else if (Ext == TEXT("cif")) 
        ParseMMCIF(CurrentPDBContent);
    else if (Ext == TEXT("sdf") || Ext == TEXT("mol"))
        ParseSDF(CurrentPDBContent);
    else 
        UE_LOG(LogTemp, Error, TEXT("Unsupported format: %s"), *Ext);
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
    FString FileTypes = TEXT("All Supported (*.pdb;*.cif;*.sdf;*.mol)|*.pdb;*.cif;*.sdf;*.mol|")
                        TEXT("PDB Files (*.pdb)|*.pdb|")
                        TEXT("mmCIF Files (*.cif)|*.cif|")
                        TEXT("SDF Files (*.sdf)|*.sdf|")
                        TEXT("MOL Files (*.mol)|*.mol");
    
    FString DefaultFile = bSave ? TEXT("structure.pdb") : TEXT("");
    FString DialogTitle = bSave ? TEXT("Save Structure") : TEXT("Load Structure");
    
    bool bOK = bSave ?
        DP->SaveFileDialog(Handle, DialogTitle, FPaths::ProjectSavedDir(), DefaultFile,
            FileTypes, EFileDialogFlags::None, Files) :
        DP->OpenFileDialog(Handle, DialogTitle, FPaths::ProjectSavedDir(), TEXT(""),
            FileTypes, EFileDialogFlags::None, Files);

    if (bOK && Files.Num() > 0) { Out = MoveTemp(Files[0]); return true; }
    return false;
}

void APDBViewer::OpenSaveDialog() { FString P; if (ShowFileDialog(true, P)) SaveStructureToFile(P); }
void APDBViewer::OpenLoadDialog() { FString P; if (ShowFileDialog(false, P)) LoadStructureFromFile(P); }

// ListView Functions for SDF Molecules

TArray<UPDBMoleculeNode*> APDBViewer::GetMoleculeNodes()
{
    TArray<UPDBMoleculeNode*> Nodes;
    TArray<FString> Keys;
    LigandMap.GetKeys(Keys);
    
    // Sort by molecule name
    Keys.Sort();
    
    for (const FString& Key : Keys)
    {
        const auto* InfoPtr = LigandMap.Find(Key);
        if (!InfoPtr || !*InfoPtr) continue;
        
        const FLigandInfo* Info = *InfoPtr;
        
        UPDBMoleculeNode* Node = NewObject<UPDBMoleculeNode>(this);
        Node->Initialize(
            Info->LigandName,
            Key,
            Info->bIsVisible,
            Info->AtomMeshes.Num(),
            Info->BondMeshes.Num()
        );
        
        Nodes.Add(Node);
    }
    
    return Nodes;
}

void APDBViewer::PopulateMoleculeListView(UListView* ListView)
{
    if (!ListView) return;
    
    // Clear existing items
    ListView->ClearListItems();
    
    // Get all molecule nodes
    TArray<UPDBMoleculeNode*> MoleculeNodes = GetMoleculeNodes();
    
    // Add each molecule to the list view
    for (UPDBMoleculeNode* MoleculeNode : MoleculeNodes)
    {
        ListView->AddItem(MoleculeNode);
    }
    
    // Request refresh
    ListView->RegenerateAllEntries();
}

void APDBViewer::ToggleMoleculeVisibility(const FString& MoleculeKey)
{
    auto** InfoPtr = LigandMap.Find(MoleculeKey);
    if (!InfoPtr || !*InfoPtr) return;
    
    auto* Info = *InfoPtr;
    Info->bIsVisible = !Info->bIsVisible;
    
    for (auto* M : Info->AtomMeshes) 
        if (M) M->SetVisibility(Info->bIsVisible);
    
    for (auto* M : Info->BondMeshes) 
        if (M) M->SetVisibility(Info->bIsVisible);
}

void APDBViewer::ToggleMoleculeNodeVisibility(UPDBMoleculeNode* Node)
{
    if (!Node) return;
    
    ToggleMoleculeVisibility(Node->MoleculeKey);
    
    // Update the node's visibility state
    auto** InfoPtr = LigandMap.Find(Node->MoleculeKey);
    if (InfoPtr && *InfoPtr)
    {
        Node->bIsVisible = (*InfoPtr)->bIsVisible;
    }
}