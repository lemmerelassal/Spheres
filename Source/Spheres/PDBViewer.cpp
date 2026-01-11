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
#include "Kismet/GameplayStatics.h"
#include "MDControlWidget.h"
#include "HydrogenGenerator.h"

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

// In PDBViewer.cpp, modify BeginPlay:

void APDBViewer::BeginPlay()
{
    Super::BeginPlay();
    FetchAndDisplayStructure(TEXT("5ENB"));

    if (auto *Cam = GetWorld()->SpawnActor<APDBCameraComponent>(APDBCameraComponent::StaticClass(), GetActorLocation(), FRotator::ZeroRotator))
        Cam->SetTargetActor(this);

    // Spawn MD control widget if class is provided
    if (MDControlWidgetClass)
    {
        APlayerController* PC = UGameplayStatics::GetPlayerController(GetWorld(), 0);
        if (PC)
        {
            MDControlWidgetInstance = CreateWidget<UMDControlWidget>(PC, MDControlWidgetClass);
            if (MDControlWidgetInstance)
            {
                MDControlWidgetInstance->AddToViewport();
                UE_LOG(LogTemp, Log, TEXT("PDBViewer: MD control widget added to viewport"));
            }
        }
    }
    
    // ADD THIS: Delay to let structure load, then test hydrogens
    FTimerHandle TimerHandle;
    GetWorld()->GetTimerManager().SetTimer(TimerHandle, [this]()
    {
        UE_LOG(LogTemp, Warning, TEXT("===== TESTING HYDROGEN GENERATION ====="));
        DebugPrintLigandInfo();
        AddExplicitHydrogens();
        int32 Count = GetHydrogenCount();
        UE_LOG(LogTemp, Warning, TEXT("Created %d hydrogen atoms"), Count);
    }, 3.0f, false);  // Wait 3 seconds for structure to load
}

void APDBViewer::FetchAndDisplayStructure(const FString &ID)
{
    CurrentStructureID = ID;
    FString URL = FString::Printf(TEXT("https://files.rcsb.org/download/%s.pdb"), *ID);
    FetchFileAsync(URL, [this, ID](bool bOK, const FString &Content)
                   {
        if (bOK) ParsePDB(Content);
        else FetchFileAsync(FString::Printf(TEXT("https://files.rcsb.org/download/%s.cif"), *ID),
            [this](bool bOK2, const FString& C) { if (bOK2) ParseMMCIF(C); }); });
}

void APDBViewer::FetchFileAsync(const FString &URL, TFunction<void(bool, const FString &)> CB)
{
    auto Req = FHttpModule::Get().CreateRequest();
    Req->SetURL(URL);
    Req->SetVerb(TEXT("GET"));
    Req->OnProcessRequestComplete().BindLambda([CB](FHttpRequestPtr R, FHttpResponsePtr Resp, bool bOK)
                                               { CB(bOK && Resp.IsValid() && Resp->GetResponseCode() == 200, bOK ? Resp->GetContentAsString() : TEXT("")); });
    Req->ProcessRequest();
}

void APDBViewer::ParsePDB(const FString &Content)
{
    CurrentPDBContent = Content;
    ClearResidueMap();
    ChainIDs.Empty();

    TArray<FString> Lines;
    Content.ParseIntoArrayLines(Lines);

    TMap<FString, TMap<FString, FVector>> ResAtoms;
    TMap<FString, FResidueMetadata> ResMeta;

    for (const auto &L : Lines)
    {
        FString RecordType = L.Mid(0, 6);
        if (L.Len() < 80 || !(RecordType.StartsWith("ATOM") || RecordType.StartsWith("HETATM")))
            continue;

        FString Chain = L.Mid(21, 1).TrimStartAndEnd();
        if (Chain.IsEmpty())
            Chain = TEXT("_");

        ChainIDs.Add(Chain);

        FString Key = FString::Printf(TEXT("%s_%s_%s"),
                                      *L.Mid(17, 3).TrimStartAndEnd(),
                                      *L.Mid(22, 4).TrimStartAndEnd(),
                                      *Chain);

        ResAtoms.FindOrAdd(Key).Add(L.Mid(12, 4).TrimStartAndEnd(),
                                    FVector(FCString::Atof(*L.Mid(30, 8)),
                                            FCString::Atof(*L.Mid(38, 8)),
                                            FCString::Atof(*L.Mid(46, 8))) *
                                        PDB::SCALE);

        if (!ResMeta.Contains(Key))
        {
            auto &M = ResMeta.Add(Key);
            M.ResidueName = L.Mid(17, 3).TrimStartAndEnd();
            M.ResidueSeq = L.Mid(22, 4).TrimStartAndEnd();
            M.Chain = Chain;
            M.RecordType = RecordType.TrimStartAndEnd();
        }
    }

    CreateResiduesFromAtomData(ResAtoms, ResMeta);
    OnResiduesLoaded.Broadcast();
    OnLigandsLoaded.Broadcast();
}

void APDBViewer::ParseMMCIF(const FString &Content)
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

    for (const auto &L : Lines)
    {
        if (L.StartsWith(TEXT("loop_")))
        {
            bLoop = true;
            Hdrs.Empty();
            continue;
        }
        if (bLoop && L.StartsWith(TEXT("_atom_site.")))
        {
            int32 I = Hdrs.Add(L);
            if (L.Contains(TEXT("Cartn_x")))
                XI = I;
            else if (L.Contains(TEXT("Cartn_y")))
                YI = I;
            else if (L.Contains(TEXT("Cartn_z")))
                ZI = I;
            else if (L.Contains(TEXT("label_comp_id")))
                RI = I;
            else if (L.Contains(TEXT("label_atom_id")))
                AI = I;
            else if (L.Contains(TEXT("group_PDB")))
                GI = I;
            else if (L.Contains(TEXT("label_asym_id")))
                CI = I;
            else if (L.Contains(TEXT("label_seq_id")) || L.Contains(TEXT("auth_seq_id")))
                SI = I;
            continue;
        }
        if (bLoop && !L.StartsWith(TEXT("_")))
        {
            TArray<FString> T;
            L.ParseIntoArrayWS(T);
            if (T.Num() > FMath::Max3(XI, YI, ZI))
                AtomTab.Add(MoveTemp(T));
        }
    }

    TMap<FString, TMap<FString, FVector>> ResAtoms;
    TMap<FString, FResidueMetadata> ResMeta;

    for (const auto &R : AtomTab)
    {
        if (XI < 0 || YI < 0 || ZI < 0 || RI < 0 || AI < 0)
            continue;

        FString GroupPDB = (GI >= 0 && R.IsValidIndex(GI)) ? R[GI] : TEXT("ATOM");
        FString Chain = (CI >= 0 && R.IsValidIndex(CI)) ? R[CI] : TEXT("_");
        if (Chain.IsEmpty())
            Chain = TEXT("_");

        ChainIDs.Add(Chain);

        FString Seq = (SI >= 0 && R.IsValidIndex(SI)) ? R[SI] : TEXT("0");
        FString Key = FString::Printf(TEXT("%s_%s_%s"), *R[RI], *Seq, *Chain);

        ResAtoms.FindOrAdd(Key).Add(R[AI],
                                    FVector(FCString::Atof(*R[XI]),
                                            FCString::Atof(*R[YI]),
                                            FCString::Atof(*R[ZI])) *
                                        PDB::SCALE);

        if (!ResMeta.Contains(Key))
        {
            auto &M = ResMeta.Add(Key);
            M.ResidueName = R[RI];
            M.ResidueSeq = Seq;
            M.Chain = Chain;
            M.RecordType = GroupPDB;
        }
    }

    CreateResiduesFromAtomData(ResAtoms, ResMeta);
    OnResiduesLoaded.Broadcast();
    OnLigandsLoaded.Broadcast();
}

void APDBViewer::ParseSDF(const FString &Content)
{
    CurrentPDBContent = Content;
    ClearLigandMap();

    TArray<FString> Lines;
    Content.ParseIntoArrayLines(Lines);

    int32 MoleculeIndex = 0;
    int32 LineIndex = 0;

    while (LineIndex < Lines.Num())
    {
        if (LineIndex + 3 >= Lines.Num())
            break;

        FString MoleculeName = Lines[LineIndex].TrimStartAndEnd();
        if (MoleculeName.IsEmpty())
            MoleculeName = FString::Printf(TEXT("MOL%d"), MoleculeIndex + 1);

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
            LineIndex++;
            continue;
        }

        FString CountsLine = Lines[CountsLineIndex];

        if (CountsLine.Len() < 6)
        {
            LineIndex++;
            continue;
        }

        int32 NumAtoms = FCString::Atoi(*CountsLine.Mid(0, 3).TrimStartAndEnd());
        int32 NumBonds = FCString::Atoi(*CountsLine.Mid(3, 3).TrimStartAndEnd());

        TArray<FVector> AtomPositions;
        TArray<FString> AtomElements;
        int32 AtomStartLine = CountsLineIndex + 1;

        for (int32 i = 0; i < NumAtoms; ++i)
        {
            int32 CurrentLine = AtomStartLine + i;
            if (CurrentLine >= Lines.Num())
                break;

            FString Line = Lines[CurrentLine];
            if (Line.Len() < 34)
                continue;

            float X = FCString::Atof(*Line.Mid(0, 10).TrimStartAndEnd());
            float Y = FCString::Atof(*Line.Mid(10, 10).TrimStartAndEnd());
            float Z = FCString::Atof(*Line.Mid(20, 10).TrimStartAndEnd());

            FString Element = Line.Mid(31, 3).TrimStartAndEnd();
            // Element symbols can be 1-3 characters: C, Cl, Uup
            // First character is always uppercase, rest are lowercase if present
            if (Element.Len() > 0)
            {
                // Standard approach: keep first char uppercase, and consecutive lowercase chars
                // This handles: C (1), Cl (2), Uup (3), etc.
                FString Normalized;
                Normalized += FChar::ToUpper(Element[0]);
                for (int32 j = 1; j < Element.Len(); ++j)
                {
                    if (FChar::IsAlpha(Element[j]))
                        Normalized += FChar::ToLower(Element[j]);
                }
                Element = Normalized;
            }
            FVector Pos(X, Y, Z);
            Pos *= PDB::SCALE;

            AtomPositions.Add(Pos);
            AtomElements.Add(Element);
        }

        int32 BondStartLine = AtomStartLine + NumAtoms;
        TArray<TPair<int32, int32>> BondPairs;
        TArray<int32> BondOrders;

        for (int32 i = 0; i < NumBonds && (BondStartLine + i) < Lines.Num(); ++i)
        {
            FString Line = Lines[BondStartLine + i];
            if (Line.Len() < 9)
                continue;

            int32 Atom1 = FCString::Atoi(*Line.Mid(0, 3).TrimStartAndEnd()) - 1;
            int32 Atom2 = FCString::Atoi(*Line.Mid(3, 3).TrimStartAndEnd()) - 1;
            int32 BondType = FCString::Atoi(*Line.Mid(6, 3).TrimStartAndEnd());

            if (Atom1 >= 0 && Atom1 < AtomPositions.Num() &&
                Atom2 >= 0 && Atom2 < AtomPositions.Num())
            {
                BondPairs.Add(TPair<int32, int32>(Atom1, Atom2));
                BondOrders.Add(BondType);
            }
        }

        
FString Key = FString::Printf(TEXT("%s"), *MoleculeName);

        auto *Info = new FLigandInfo();
        Info->LigandName = MoleculeName;
        Info->bIsVisible = false;
        
        // Store atom positions and elements for hydrogen generation
        Info->AtomPositions = AtomPositions;
        Info->AtomElements = AtomElements;
        
        // Store bond connectivity
        for (int32 i = 0; i < BondPairs.Num(); ++i)
        {
            Info->BondPairs.Add(BondPairs[i]);
            Info->BondOrders.Add(BondOrders[i]); // Store the bond order too
        }

        // Generate and add hydrogens immediately if auto-generation is enabled
        if (bAutoGenerateHydrogens)
        {
            TArray<TPair<FVector, int32>> Hydrogens = FHydrogenGenerator::GenerateHydrogens(
                Info->AtomPositions,
                Info->AtomElements,
                Info->BondPairs,
                Info->BondOrders);
            
            UE_LOG(LogTemp, Log, TEXT("Generated %d hydrogens for molecule %s"), Hydrogens.Num(), *MoleculeName);
            
            // Add hydrogen positions and elements to the ligand info
            for (const auto& HPair : Hydrogens)
            {
                FVector HPos = HPair.Key;
                int32 ParentIdx = HPair.Value;
                
                // Add to positions and elements arrays
                Info->AtomPositions.Add(HPos);
                Info->AtomElements.Add(TEXT("H"));
                
                // Draw hydrogen sphere
                DrawSphere(HPos.X, HPos.Y, HPos.Z, FLinearColor::White, GetRootComponent(), Info->AtomMeshes);
                
                // Draw bond to parent
                if (Info->AtomPositions.IsValidIndex(ParentIdx))
                {
                    FVector ParentPos = Info->AtomPositions[ParentIdx];
                    FString ParentElement = Info->AtomElements.IsValidIndex(ParentIdx) 
                        ? Info->AtomElements[ParentIdx] 
                        : TEXT("C");
                    
                    DrawBond(ParentPos, HPos, 1, ParentElement, TEXT("H"), GetRootComponent(), Info->BondMeshes);
                }
            }
        }

        // Draw heavy atoms
        for (int32 i = 0; i < AtomPositions.Num(); ++i)
        {
            FLinearColor Color = GetElementColor(AtomElements[i]);
            DrawSphere(AtomPositions[i].X, AtomPositions[i].Y, AtomPositions[i].Z,
                       Color, GetRootComponent(), Info->AtomMeshes);
        }

        // Draw bonds between heavy atoms
        for (int32 i = 0; i < BondPairs.Num(); ++i)
        {
            int32 Atom1 = BondPairs[i].Key;
            int32 Atom2 = BondPairs[i].Value;
            int32 Order = BondOrders[i];

            if (Order == 4)
                Order = 1;

            DrawBond(AtomPositions[Atom1], AtomPositions[Atom2], Order,
                     AtomElements[Atom1], AtomElements[Atom2], GetRootComponent(), Info->BondMeshes);
        }

        for (auto *Mesh : Info->AtomMeshes)
        {
            if (Mesh)
                Mesh->SetVisibility(Info->bIsVisible);
        }

        for (auto *Mesh : Info->BondMeshes)
        {
            if (Mesh)
                Mesh->SetVisibility(Info->bIsVisible);
        }

        LigandMap.Add(Key, Info);

        LineIndex = BondStartLine + NumBonds;
        while (LineIndex < Lines.Num())
        {
            if (Lines[LineIndex].StartsWith(TEXT("$$")))
            {
                LineIndex++;
                break;
            }
            LineIndex++;
        }

        MoleculeIndex++;
    }

    OnLigandsLoaded.Broadcast();
}

void APDBViewer::CreateResiduesFromAtomData(const TMap<FString, TMap<FString, FVector>> &ResAtoms, const TMap<FString, FResidueMetadata> &Meta)
{
    for (const auto &P : ResAtoms)
    {
        const auto *M = Meta.Find(P.Key);
        if (!M)
            continue;

        // Check if this is a HETATM entry - add to ligand map instead
        
        
        if (M->RecordType.StartsWith(TEXT("HETATM")))
        {
            auto *LigInfo = new FLigandInfo();
            LigInfo->LigandName = FString::Printf(TEXT("%s-%s-%s"), *M->ResidueName, *M->ResidueSeq, *M->Chain);
            LigInfo->bIsVisible = false;

            // Store atom positions and elements
            for (const auto &A : P.Value)
            {
                FString Element = A.Key.TrimStartAndEnd().Left(1);
                if (A.Key.Len() > 1 && FChar::IsLower(A.Key[1]))
                    Element = A.Key.Left(2);
                
                LigInfo->AtomPositions.Add(A.Value);
                LigInfo->AtomElements.Add(Element);
                
                DrawSphere(A.Value.X, A.Value.Y, A.Value.Z, GetElementColor(Element), GetRootComponent(), LigInfo->AtomMeshes);
            }

            // Set initial visibility
            for (auto *Mesh : LigInfo->AtomMeshes)
            {
                if (Mesh)
                    Mesh->SetVisibility(LigInfo->bIsVisible);
            }

            LigandMap.Add(P.Key, LigInfo);
            FetchLigandBondsForHETATM(P.Key, M->ResidueName, P.Value, LigInfo);
        }
        else
        {
            // Regular ATOM entry - add to residue map
            auto *Info = new FResidueInfo();
            Info->ResidueName = M->ResidueName;
            Info->ResidueSeq = M->ResidueSeq;
            Info->Chain = M->Chain;
            Info->RecordType = M->RecordType;
            Info->bIsVisible = true;

            // Store raw atom positions/elements for MMGBSA and other analysis
            for (const auto &A : P.Value)
            {
                FString Element = A.Key.TrimStartAndEnd().Left(1);
                if (A.Key.Len() > 1 && FChar::IsLower(A.Key[1]))
                    Element = A.Key.Left(2);
                Info->AtomPositions.Add(A.Value);
                Info->AtomElements.Add(Element);
            }

            // For ATOM records, only draw bonds (no spheres)
            DrawProteinBonds(P.Value, Info);

            ResidueMap.Add(P.Key, Info);
        }
    }
}

void APDBViewer::DrawProteinBonds(const TMap<FString, FVector> &AtomPositions, FResidueInfo *ResInfo)
{
    if (!ResInfo)
        return;

    const float BondThreshold = 100.0f;

    TArray<TPair<FString, FVector>> Atoms;
    for (const auto &Pair : AtomPositions)
    {
        Atoms.Add(TPair<FString, FVector>(Pair.Key, Pair.Value));
    }

    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Distance = FVector::Dist(Atoms[i].Value, Atoms[j].Value);

            if (Distance < BondThreshold)
            {
                FString Element1 = Atoms[i].Key.TrimStartAndEnd().Left(1);
                if (Atoms[i].Key.Len() > 1 && FChar::IsLower(Atoms[i].Key[1]))
                    Element1 = Atoms[i].Key.Left(2);

                FString Element2 = Atoms[j].Key.TrimStartAndEnd().Left(1);
                if (Atoms[j].Key.Len() > 1 && FChar::IsLower(Atoms[j].Key[1]))
                    Element2 = Atoms[j].Key.Left(2);

                DrawBond(Atoms[i].Value, Atoms[j].Value, 1, Element1, Element2, GetRootComponent(), ResInfo->BondMeshes);
            }
        }
    }
}

void APDBViewer::FetchLigandBondsForHETATM(const FString &Key, const FString &Name, const TMap<FString, FVector> &Pos, FLigandInfo *LigInfo)
{
    FetchFileAsync(FString::Printf(TEXT("https://files.rcsb.org/ligands/download/%s.cif"), *Name.ToUpper()),
                   [this, Key, Pos, LigInfo](bool bOK, const FString &C)
                   {
                       if (bOK && LigandMap.Contains(Key))
                           ParseLigandCIFForLigand(C, Pos, LigInfo);
                   });
}


void APDBViewer::ParseLigandCIFForLigand(const FString &Content, const TMap<FString, FVector> &Pos, FLigandInfo *Info)
{
    if (!Info)
        return;

    TMap<FString, FVector> NormPos;
    TMap<FString, FString> AtomElements;

    for (const auto &P : Pos)
    {
        FString K;
        for (const TCHAR C : P.Key)
            if (FChar::IsAlnum(C))
                K.AppendChar(FChar::ToUpper(C));
        if (!K.IsEmpty())
        {
            NormPos.Add(K, P.Value);
            FString Element = P.Key.TrimStartAndEnd().Left(1);
            if (P.Key.Len() > 1 && FChar::IsLower(P.Key[1]))
                Element = P.Key.Left(2);
            AtomElements.Add(K, Element);
        }
    }

    TArray<FString> Lines;
    Content.ParseIntoArrayLines(Lines);

    TArray<FString> Hdrs;
    int32 A1 = -1, A2 = -1, BO = -1;
    bool bLoop = false;

    for (const auto &L : Lines)
    {
        if (L.StartsWith(TEXT("loop_")))
        {
            bLoop = true;
            Hdrs.Empty();
            A1 = A2 = BO = -1;
            continue;
        }
        if (bLoop && L.StartsWith(TEXT("_")))
        {
            int32 I = Hdrs.Add(L);
            FString Lo = L.ToLower();
            if (Lo.Contains(TEXT("atom_id_1")) || Lo.Contains(TEXT("atom_1")))
                A1 = I;
            else if (Lo.Contains(TEXT("atom_id_2")) || Lo.Contains(TEXT("atom_2")))
                A2 = I;
            else if (Lo.Contains(TEXT("value_order")) || Lo.Contains(TEXT("bond_order")))
                BO = I;
            continue;
        }
        if (bLoop && !L.StartsWith(TEXT("_")) && !L.StartsWith(TEXT("data_")))
        {
            if (A1 < 0 || A2 < 0)
                continue;
            TArray<FString> T;
            L.ParseIntoArrayWS(T);
            if (T.Num() <= FMath::Max(A1, A2))
                continue;

            FString ID1 = NormalizeAtomID(T[A1]);
            FString ID2 = NormalizeAtomID(T[A2]);
            int32 Ord = ParseBondOrder(BO >= 0 && T.IsValidIndex(BO) ? T[BO] : TEXT("1"));

            const auto *P1 = NormPos.Find(ID1);
            const auto *P2 = NormPos.Find(ID2);
            if (P1 && P2)
            {
                // Find atom indices in the stored positions
                int32 Idx1 = -1, Idx2 = -1;
                for (int32 i = 0; i < Info->AtomPositions.Num(); ++i)
                {
                    if (Info->AtomPositions[i].Equals(*P1, 0.1f))
                        Idx1 = i;
                    if (Info->AtomPositions[i].Equals(*P2, 0.1f))
                        Idx2 = i;
                }
                
                // Store bond connectivity
                if (Idx1 >= 0 && Idx2 >= 0)
                {
                    Info->BondPairs.Add(TPair<int32, int32>(Idx1, Idx2));
                    Info->BondOrders.Add(Ord);
                }
                
                // Extract element symbols from atom IDs (first 1-2 characters)
                FString Elem1, Elem2;
                if (ID1.Len() > 0)
                {
                    Elem1 = ID1.Left(1);
                    if (ID1.Len() > 1 && FChar::IsLower(ID1[1]))
                        Elem1 = ID1.Left(2);
                }
                if (ID2.Len() > 0)
                {
                    Elem2 = ID2.Left(1);
                    if (ID2.Len() > 1 && FChar::IsLower(ID2[1]))
                        Elem2 = ID2.Left(2);
                }
                
                DrawBond(*P1, *P2, Ord, Elem1, Elem2, GetRootComponent(), Info->BondMeshes);

                if (Info->BondMeshes.Num() > 0)
                {
                    Info->BondMeshes.Last()->SetVisibility(Info->bIsVisible);
                }
            }
        }
        else if (bLoop && L.StartsWith(TEXT("data_")))
            break;
    }
}

FString APDBViewer::NormalizeAtomID(const FString &In) const
{
    FString Out;
    for (const TCHAR C : In)
        if (FChar::IsAlnum(C))
            Out.AppendChar(FChar::ToUpper(C));
    return Out;
}

int32 APDBViewer::ParseBondOrder(const FString &S) const
{
    if (S.Len() == 1 && S[0] >= '1' && S[0] <= '3')
        return S[0] - '0';
    FString L = S.ToLower();
    if (L.Contains(TEXT("ar")))
        return 1;
    if (L.Contains(TEXT("doub")) || L.Equals(TEXT("d")))
        return 2;
    if (L.Contains(TEXT("trip")) || L.Equals(TEXT("t")))
        return 3;
    return FMath::Clamp(FCString::Atoi(*S), 1, 3);
}

void APDBViewer::DrawSphere(float X, float Y, float Z, const FLinearColor &Col, USceneComponent *Par, TArray<UStaticMeshComponent *> &Out)
{
    if (!SphereMeshAsset || !SphereMaterialAsset || !Par)
        return;

    auto *Sph = NewObject<UStaticMeshComponent>(this);
    Sph->SetStaticMesh(SphereMeshAsset);
    Sph->SetWorldScale3D(FVector(PDB::SPHERE_SIZE));
    Sph->SetCollisionEnabled(ECollisionEnabled::NoCollision);

    auto *Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, Sph);
    Mat->SetVectorParameterValue(TEXT("Color"), Col);
    Mat->SetScalarParameterValue(TEXT("EmissiveIntensity"), 5.0f);
    Sph->SetMaterial(0, Mat);

    Sph->AttachToComponent(Par, FAttachmentTransformRules::KeepWorldTransform);
    Sph->SetRelativeLocation(Par->GetComponentTransform().InverseTransformPosition(FVector(X, Y, Z)));
    Sph->RegisterComponent();

    Out.Add(Sph);
    AllAtomMeshes.Add(Sph);
}

void APDBViewer::DrawBond(const FVector &S, const FVector &E, int32 Ord, const FString &Element1, const FString &Element2, USceneComponent *Par, TArray<UStaticMeshComponent *> &Out)
{
    if (!CylinderMeshAsset || !SphereMaterialAsset || !Par)
        return;

    FVector V = E - S;
    float Len = V.Size();
    if (Len < KINDA_SMALL_NUMBER)
        return;

    FRotator Rot = FRotationMatrix::MakeFromZ(V).Rotator();
    float Scale = Len / 100.0f;

    FLinearColor Color1 = GetElementColor(Element1);
    FLinearColor Color2 = GetElementColor(Element2);

    auto MakeCyl = [&](const FVector &Pos, float ScaleZ, const FLinearColor &Color)
    {
        auto *Cyl = NewObject<UStaticMeshComponent>(this);
        Cyl->SetStaticMesh(CylinderMeshAsset);
        Cyl->SetWorldLocation(Pos);
        Cyl->SetWorldRotation(Rot);
        Cyl->SetWorldScale3D(FVector(PDB::CYLINDER_SIZE, PDB::CYLINDER_SIZE, ScaleZ));
        Cyl->SetCollisionEnabled(ECollisionEnabled::NoCollision);
        auto *Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, Cyl);
        Mat->SetVectorParameterValue(TEXT("Color"), Color);
        Mat->SetScalarParameterValue(TEXT("EmissiveIntensity"), 2.0f);
        Cyl->SetMaterial(0, Mat);
        Cyl->AttachToComponent(Par, FAttachmentTransformRules::KeepWorldTransform);
        Cyl->RegisterComponent();
        Out.Add(Cyl);
        AllBondMeshes.Add(Cyl);
    };

    FVector HalfDir = V * 0.5f;
    FVector FirstMid = S + V * 0.25f;
    FVector SecondMid = S + V * 0.75f;
    float HalfScale = Scale * 0.5f;

    if (Ord <= 1)
    {
        MakeCyl(FirstMid, HalfScale, Color1);
        MakeCyl(SecondMid, HalfScale, Color2);
        return;
    }

    FVector Perp = FVector::CrossProduct(V.GetSafeNormal(), FVector::UpVector);
    if (Perp.SizeSquared() < KINDA_SMALL_NUMBER)
        Perp = FVector::CrossProduct(V.GetSafeNormal(), FVector::RightVector);
    Perp.Normalize();
    FVector Off = Perp * PDB::BOND_OFFSET;

    if (Ord == 2)
    {
        MakeCyl(FirstMid + Off, HalfScale, Color1);
        MakeCyl(SecondMid + Off, HalfScale, Color2);
        MakeCyl(FirstMid - Off, HalfScale, Color1);
        MakeCyl(SecondMid - Off, HalfScale, Color2);
    }
    else if (Ord == 3)
    {
        MakeCyl(FirstMid, HalfScale, Color1);
        MakeCyl(SecondMid, HalfScale, Color2);
        MakeCyl(FirstMid + Off, HalfScale, Color1);
        MakeCyl(SecondMid + Off, HalfScale, Color2);
        MakeCyl(FirstMid - Off, HalfScale, Color1);
        MakeCyl(SecondMid - Off, HalfScale, Color2);
    }
    else
    {
        MakeCyl(FirstMid, HalfScale, Color1);
        MakeCyl(SecondMid, HalfScale, Color2);
    }
}


void APDBViewer::ClearOverlapMarkers()
{
    for (auto *M : OverlapMarkers)
    {
        if (M && IsValid(M))
            M->DestroyComponent();
    }
    OverlapMarkers.Empty();
}


void APDBViewer::HighlightOverlapAtoms(const TArray<FMMOverlapInfo>& Overlaps)
{
    ClearOverlapMarkers();
    if (!SphereMeshAsset || !SphereMaterialAsset || !CylinderMeshAsset)
        return;

    APlayerController* PC = UGameplayStatics::GetPlayerController(GetWorld(), 0);
    FVector FirstPos = FVector::ZeroVector;
    bool bHaveFirst = false;

    for (const FMMOverlapInfo& O : Overlaps)
    {
        // Create small red sphere at receptor position
        auto *S1 = NewObject<UStaticMeshComponent>(this);
        S1->SetStaticMesh(SphereMeshAsset);
        S1->SetWorldScale3D(FVector(PDB::SPHERE_SIZE * 0.6f));
        S1->SetCollisionEnabled(ECollisionEnabled::NoCollision);
        auto *Mat1 = UMaterialInstanceDynamic::Create(SphereMaterialAsset, S1);
        Mat1->SetVectorParameterValue(TEXT("Color"), FLinearColor::Red);
        Mat1->SetScalarParameterValue(TEXT("EmissiveIntensity"), 8.0f);
        S1->SetMaterial(0, Mat1);
        S1->AttachToComponent(GetRootComponent(), FAttachmentTransformRules::KeepWorldTransform);
        S1->SetWorldLocation(O.ReceptorPosition);
        S1->RegisterComponent();
        OverlapMarkers.Add(S1);

        // Create small red sphere at ligand position
        auto *S2 = NewObject<UStaticMeshComponent>(this);
        S2->SetStaticMesh(SphereMeshAsset);
        S2->SetWorldScale3D(FVector(PDB::SPHERE_SIZE * 0.6f));
        S2->SetCollisionEnabled(ECollisionEnabled::NoCollision);
        auto *Mat2 = UMaterialInstanceDynamic::Create(SphereMaterialAsset, S2);
        Mat2->SetVectorParameterValue(TEXT("Color"), FLinearColor::Red);
        Mat2->SetScalarParameterValue(TEXT("EmissiveIntensity"), 8.0f);
        S2->SetMaterial(0, Mat2);
        S2->AttachToComponent(GetRootComponent(), FAttachmentTransformRules::KeepWorldTransform);
        S2->SetWorldLocation(O.LigandPosition);
        S2->RegisterComponent();
        OverlapMarkers.Add(S2);

        // Create a thin cylinder connecting the two
        FVector Dir = O.LigandPosition - O.ReceptorPosition;
        float Len = Dir.Size();
        if (Len > KINDA_SMALL_NUMBER)
        {
            FRotator Rot = FRotationMatrix::MakeFromZ(Dir).Rotator();
            FVector Mid = O.ReceptorPosition + 0.5f * Dir;
            auto *Cyl = NewObject<UStaticMeshComponent>(this);
            Cyl->SetStaticMesh(CylinderMeshAsset);
            Cyl->SetWorldLocation(Mid);
            Cyl->SetWorldRotation(Rot);
            Cyl->SetWorldScale3D(FVector(PDB::CYLINDER_SIZE * 0.4f, PDB::CYLINDER_SIZE * 0.4f, Len / 100.0f));
            Cyl->SetCollisionEnabled(ECollisionEnabled::NoCollision);
            auto *MatC = UMaterialInstanceDynamic::Create(SphereMaterialAsset, Cyl);
            MatC->SetVectorParameterValue(TEXT("Color"), FLinearColor::Red);
            MatC->SetScalarParameterValue(TEXT("EmissiveIntensity"), 6.0f);
            Cyl->SetMaterial(0, MatC);
            Cyl->AttachToComponent(GetRootComponent(), FAttachmentTransformRules::KeepWorldTransform);
            Cyl->RegisterComponent();
            OverlapMarkers.Add(Cyl);
        }

        if (!bHaveFirst)
        {
            FirstPos = O.LigandPosition;
            bHaveFirst = true;
        }
    }

    // Center camera on first overlap (optional)
    if (bHaveFirst && PC)
    {
        PC->SetViewTarget(this);
    }
}

FLinearColor APDBViewer::GetElementColor(const FString &E) const
{
    static const TMap<FString, FLinearColor> Colors = {
        {TEXT("C"), FLinearColor(0.1f, 0.1f, 0.1f)}, {TEXT("O"), FLinearColor::Red}, {TEXT("H"), FLinearColor::White}, {TEXT("D"), FLinearColor::White}, {TEXT("N"), FLinearColor::Blue}, {TEXT("S"), FLinearColor::Yellow}, {TEXT("CL"), FLinearColor(0, 1, 0)}, {TEXT("P"), FLinearColor(1, 0.5f, 0)}, {TEXT("F"), FLinearColor(0, 1, 0)}, {TEXT("BR"), FLinearColor(0.6f, 0.2f, 0.2f)}, {TEXT("I"), FLinearColor(0.4f, 0, 0.8f)}, {TEXT("FE"), FLinearColor(0.8f, 0.4f, 0)}, {TEXT("MG"), FLinearColor(0, 0.8f, 0)}, {TEXT("ZN"), FLinearColor(0.5f, 0.5f, 0.5f)}, {TEXT("CA"), FLinearColor(0.2f, 0.6f, 1)}, {TEXT("NA"), FLinearColor(0, 0, 1)}, {TEXT("K"), FLinearColor(0.5f, 0, 1)}, {TEXT("CU"), FLinearColor(1, 0.5f, 0)}, {TEXT("B"), FLinearColor(1, 0.7f, 0.7f)}};
    const auto *C = Colors.Find(E.ToUpper());
    return C ? *C : FLinearColor::Gray;
}

void APDBViewer::ClearResidueMap()
{
    for (auto &P : ResidueMap)
        delete P.Value;
    ResidueMap.Empty();
    ChainIDs.Empty();
}

void APDBViewer::ClearLigandMap()
{
    for (auto &P : LigandMap)
        delete P.Value;
    LigandMap.Empty();
}


void APDBViewer::ClearCurrentStructure()
{
    ClearResidueMap();
    ClearHydrogens(); // Add this line
    
    for (auto *M : AllAtomMeshes)
        if (M && IsValid(M))
            M->DestroyComponent();
    for (auto *M : AllBondMeshes)
        if (M && IsValid(M))
            M->DestroyComponent();
    AllAtomMeshes.Empty();
    AllBondMeshes.Empty();
}

void APDBViewer::ToggleResidueVisibility(const FString &Key)
{
    auto **InfoPtr = ResidueMap.Find(Key);
    if (!InfoPtr || !*InfoPtr)
        return;
    auto *Info = *InfoPtr;
    Info->bIsVisible = !Info->bIsVisible;
    for (auto *M : Info->AtomMeshes)
        if (M)
            M->SetVisibility(Info->bIsVisible);
    for (auto *M : Info->BondMeshes)
        if (M)
            M->SetVisibility(Info->bIsVisible);
}

TArray<UPDBTreeNode *> APDBViewer::GetChainNodes()
{
    TArray<UPDBTreeNode *> Nodes;
    TArray<FString> SortedChains = ChainIDs.Array();
    SortedChains.Sort();

    for (const FString &ChainID : SortedChains)
    {
        FString DisplayName = ChainID == TEXT("_")
                                  ? TEXT("Chain (No ID)")
                                  : FString::Printf(TEXT("Chain %s"), *ChainID);

        UPDBTreeNode *Node = NewObject<UPDBTreeNode>(this);
        Node->Initialize(DisplayName, ChainID, true, ChainID);
        Nodes.Add(Node);
    }

    return Nodes;
}

TArray<UPDBTreeNode *> APDBViewer::GetResidueNodesForChain(const FString &ChainID)
{
    TArray<UPDBTreeNode *> Nodes;
    TArray<FString> Keys;
    ResidueMap.GetKeys(Keys);

    Keys.Sort([](const FString &A, const FString &B)
              {
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
        return A < B; });

    for (const FString &Key : Keys)
    {
        const auto *InfoPtr = ResidueMap.Find(Key);
        if (!InfoPtr || !*InfoPtr)
            continue;

        const FResidueInfo *Info = *InfoPtr;
        if (Info->Chain != ChainID)
            continue;

        FString DisplayName = FString::Printf(TEXT("%s %s"), *Info->ResidueName, *Info->ResidueSeq);

        UPDBTreeNode *Node = NewObject<UPDBTreeNode>(this);
        Node->Initialize(DisplayName, Key, false, ChainID);
        Node->bIsVisible = Info->bIsVisible;

        Nodes.Add(Node);
    }

    return Nodes;
}

void APDBViewer::ToggleChainVisibility(const FString &ChainID)
{
    bool bNewVisibility = false;
    bool bFirstResidue = true;

    for (const auto &Pair : ResidueMap)
    {
        if (Pair.Value && Pair.Value->Chain == ChainID)
        {
            if (bFirstResidue)
            {
                bNewVisibility = !Pair.Value->bIsVisible;
                bFirstResidue = false;
            }

            Pair.Value->bIsVisible = bNewVisibility;
            for (auto *M : Pair.Value->AtomMeshes)
                if (M)
                    M->SetVisibility(bNewVisibility);
            for (auto *M : Pair.Value->BondMeshes)
                if (M)
                    M->SetVisibility(bNewVisibility);
        }
    }
}

void APDBViewer::ToggleNodeVisibility(UPDBTreeNode *Node)
{
    if (!Node)
        return;

    if (Node->bIsChain)
    {
        ToggleChainVisibility(Node->ChainID);
    }
    else
    {
        ToggleResidueVisibility(Node->NodeKey);
    }
}

void APDBViewer::PopulateTreeView(UTreeView *TreeView)
{
    if (!TreeView)
        return;

    TreeView->ClearListItems();

    TArray<UPDBTreeNode *> ChainNodes = GetChainNodes();

    for (UPDBTreeNode *ChainNode : ChainNodes)
    {
        TreeView->AddItem(ChainNode);
    }

    TreeView->RequestRefresh();
}

TArray<UObject *> APDBViewer::GetChildrenForNode(UPDBTreeNode *Node)
{
    TArray<UObject *> ChildNodes;

    if (!Node)
        return ChildNodes;

    if (Node->bIsChain)
    {
        TArray<UPDBTreeNode *> Residues = GetResidueNodesForChain(Node->ChainID);

        for (UPDBTreeNode *Residue : Residues)
        {
            ChildNodes.Add(Residue);
        }
    }

    return ChildNodes;
}

TArray<FString> APDBViewer::GetResidueList() const
{
    TArray<FString> List;
    ResidueMap.GetKeys(List);
    List.Sort([](const FString &A, const FString &B)
              {
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
        return A < B; });
    return List;
}

TArray<FString> APDBViewer::GetLigandList() const { return GetResidueList(); }

FString APDBViewer::GetResidueDisplayName(const FString &Key) const
{
    const auto *P = ResidueMap.Find(Key);
    return (P && *P)
               ? FString::Printf(TEXT("%s %s"), *(*P)->ResidueName, *(*P)->ResidueSeq)
               : Key;
}

FString APDBViewer::GetLigandDisplayName(const FString &Key) const
{
    const auto *P = LigandMap.Find(Key);
    return (P && *P)
               ? FString::Printf(TEXT("%s"), *(*P)->LigandName)
               : Key;
}

FLigandInfo *APDBViewer::GetVisibleLigandInfo() const
{
    FLigandInfo *BestInfo = nullptr;
    int32 BestCount = 0;

    for (auto &Pair : LigandMap)
    {
        FLigandInfo *Info = Pair.Value;
        if (!Info || !Info->bIsVisible)
            continue;
        int32 Count = Info->AtomMeshes.Num() + Info->BondMeshes.Num();
        if (Count > BestCount)
        {
            BestCount = Count;
            BestInfo = Info;
        }
    }

    return BestInfo;
}

void APDBViewer::SaveStructureToFile(const FString &Path)
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

void APDBViewer::LoadStructureFromFile(const FString &Path)
{
    FString Content;
    if (!FFileHelper::LoadFileToString(Content, *Path))
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to load: %s"), *Path);
        return;
    }

    FString Ext = FPaths::GetExtension(Path).ToLower();

    if (Ext == TEXT("pdb") || Ext == TEXT("cif"))
    {
        ClearCurrentStructure();
        CurrentStructureID = FPaths::GetBaseFilename(Path);
    }

    CurrentPDBContent = Content;

    if (Ext == TEXT("pdb"))
        ParsePDB(CurrentPDBContent);
    else if (Ext == TEXT("cif"))
        ParseMMCIF(CurrentPDBContent);
    else if (Ext == TEXT("sdf") || Ext == TEXT("mol"))
        ParseSDF(CurrentPDBContent);
    else
        UE_LOG(LogTemp, Error, TEXT("Unsupported format: %s"), *Ext);
}

bool APDBViewer::ShowFileDialog(bool bSave, FString &Out)
{
    auto *DP = FDesktopPlatformModule::Get();
    if (!DP)
        return false;

    void *Handle = nullptr;
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

    bool bOK = bSave ? DP->SaveFileDialog(Handle, DialogTitle, FPaths::ProjectSavedDir(), DefaultFile,
                                          FileTypes, EFileDialogFlags::None, Files)
                     : DP->OpenFileDialog(Handle, DialogTitle, FPaths::ProjectSavedDir(), TEXT(""),
                                          FileTypes, EFileDialogFlags::None, Files);

    if (bOK && Files.Num() > 0)
    {
        Out = MoveTemp(Files[0]);
        return true;
    }
    return false;
}

void APDBViewer::OpenSaveDialog()
{
    FString P;
    if (ShowFileDialog(true, P))
        SaveStructureToFile(P);
}
void APDBViewer::OpenLoadDialog()
{
    FString P;
    if (ShowFileDialog(false, P))
        LoadStructureFromFile(P);
}

TArray<UPDBMoleculeNode *> APDBViewer::GetMoleculeNodes()
{
    TArray<UPDBMoleculeNode *> Nodes;
    TArray<FString> Keys;
    LigandMap.GetKeys(Keys);

    Keys.Sort();

    for (const FString &Key : Keys)
    {
        const auto *InfoPtr = LigandMap.Find(Key);
        if (!InfoPtr || !*InfoPtr)
            continue;

        const FLigandInfo *Info = *InfoPtr;

        UPDBMoleculeNode *Node = NewObject<UPDBMoleculeNode>(this);
        Node->Initialize(
            Info->LigandName,
            Key,
            Info->bIsVisible,
            Info->AtomMeshes.Num(),
            Info->BondMeshes.Num());

        Nodes.Add(Node);
    }

    return Nodes;
}

void APDBViewer::PopulateMoleculeListView(UListView *ListView)
{
    if (!ListView)
        return;

    ListView->ClearListItems();

    TArray<UPDBMoleculeNode *> MoleculeNodes = GetMoleculeNodes();

    for (UPDBMoleculeNode *MoleculeNode : MoleculeNodes)
    {
        ListView->AddItem(MoleculeNode);
    }

    ListView->RegenerateAllEntries();
}

void APDBViewer::ToggleMoleculeVisibility(const FString &MoleculeKey)
{
    auto **InfoPtr = LigandMap.Find(MoleculeKey);
    if (!InfoPtr || !*InfoPtr)
        return;

    auto *Info = *InfoPtr;
    bool bNewVisible = !Info->bIsVisible;

    if (bNewVisible)
    {
        for (auto &P : LigandMap)
        {
            if (P.Key != MoleculeKey && P.Value)
            {
                P.Value->bIsVisible = false;
                for (auto *M : P.Value->AtomMeshes)
                    if (M)
                        M->SetVisibility(false);
                for (auto *M : P.Value->BondMeshes)
                    if (M)
                        M->SetVisibility(false);
            }
        }
    }

    Info->bIsVisible = bNewVisible;

    for (auto *M : Info->AtomMeshes)
        if (M)
            M->SetVisibility(Info->bIsVisible);

    for (auto *M : Info->BondMeshes)
        if (M)
            M->SetVisibility(Info->bIsVisible);

    OnLigandsLoaded.Broadcast();
}

void APDBViewer::ToggleMoleculeNodeVisibility(UPDBMoleculeNode *Node)
{
    if (!Node)
        return;

    ToggleMoleculeVisibility(Node->MoleculeKey);

    auto **InfoPtr = LigandMap.Find(Node->MoleculeKey);
    if (InfoPtr && *InfoPtr)
    {
        Node->bIsVisible = (*InfoPtr)->bIsVisible;
    }
}


void APDBViewer::AddExplicitHydrogens()
{
    if (bHydrogensVisible)
    {
        UE_LOG(LogTemp, Warning, TEXT("Hydrogens already visible"));
        return;
    }
    
    int32 TotalHydrogens = 0;
    
    // Add hydrogens to all visible ligands
    for (auto& Pair : LigandMap)
    {
        FLigandInfo* LigInfo = Pair.Value;
        if (LigInfo && LigInfo->bIsVisible)
        {
            UE_LOG(LogTemp, Log, TEXT("Adding hydrogens to ligand: %s"), *Pair.Key);
            UE_LOG(LogTemp, Log, TEXT("  Atoms: %d, Elements: %d, Bonds: %d"), 
                LigInfo->AtomPositions.Num(), 
                LigInfo->AtomElements.Num(),
                LigInfo->BondPairs.Num());
            
            int32 Count = AddHydrogensToLigand(LigInfo);
            TotalHydrogens += Count;
            UE_LOG(LogTemp, Log, TEXT("  Generated %d hydrogens"), Count);
        }
    }
    
    bHydrogensVisible = true;
    UE_LOG(LogTemp, Log, TEXT("Added %d total explicit hydrogens"), TotalHydrogens);
}


void APDBViewer::RemoveExplicitHydrogens()
{
    if (!bHydrogensVisible)
        return;
    
    ClearHydrogens();
    bHydrogensVisible = false;
    UE_LOG(LogTemp, Log, TEXT("Removed explicit hydrogens"));
}

void APDBViewer::ToggleHydrogens()
{
    if (bHydrogensVisible)
        RemoveExplicitHydrogens();
    else
        AddExplicitHydrogens();
}

int32 APDBViewer::AddHydrogensToLigand(FLigandInfo* LigInfo)
{
    if (!LigInfo)
    {
        UE_LOG(LogTemp, Warning, TEXT("LigInfo is null"));
        return 0;
    }
    
    if (LigInfo->AtomPositions.Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("No atom positions stored in ligand"));
        return 0;
    }
    
    if (LigInfo->AtomElements.Num() != LigInfo->AtomPositions.Num())
    {
        UE_LOG(LogTemp, Warning, TEXT("Atom positions (%d) and elements (%d) count mismatch"), 
            LigInfo->AtomPositions.Num(), LigInfo->AtomElements.Num());
        return 0;
    }
    
    // Use stored bond pairs and orders if available, otherwise infer from distances
    TArray<TPair<int32, int32>> BondPairs = LigInfo->BondPairs;
    TArray<int32> BondOrders = LigInfo->BondOrders;
    
    if (BondPairs.Num() == 0)
    {
        UE_LOG(LogTemp, Log, TEXT("No bond pairs stored, inferring from distances..."));
        const float BondThreshold = 100.0f; // Scaled distance threshold
        
        for (int32 i = 0; i < LigInfo->AtomPositions.Num(); ++i)
        {
            for (int32 j = i + 1; j < LigInfo->AtomPositions.Num(); ++j)
            {
                float Dist = FVector::Dist(LigInfo->AtomPositions[i], LigInfo->AtomPositions[j]);
                if (Dist < BondThreshold)
                {
                    BondPairs.Add(TPair<int32, int32>(i, j));
                    BondOrders.Add(1); // Assume single bond when inferring
                }
            }
        }
        UE_LOG(LogTemp, Log, TEXT("Inferred %d bonds"), BondPairs.Num());
    }
    
    // Ensure bond orders array matches bond pairs
    if (BondOrders.Num() != BondPairs.Num())
    {
        UE_LOG(LogTemp, Warning, TEXT("Bond orders (%d) don't match bond pairs (%d), filling with 1s"), 
            BondOrders.Num(), BondPairs.Num());
        BondOrders.SetNum(BondPairs.Num());
        for (int32 i = 0; i < BondOrders.Num(); ++i)
        {
            if (BondOrders[i] == 0)
                BondOrders[i] = 1;
        }
    }
    
    // Generate hydrogen positions with bond orders
    TArray<TPair<FVector, int32>> Hydrogens = FHydrogenGenerator::GenerateHydrogens(
        LigInfo->AtomPositions,
        LigInfo->AtomElements,
        BondPairs,
        BondOrders);
    
    UE_LOG(LogTemp, Log, TEXT("Generated %d hydrogen positions"), Hydrogens.Num());
    
    // Draw hydrogen atoms
    for (const auto& HPair : Hydrogens)
    {
        FVector HPos = HPair.Key;
        int32 ParentIdx = HPair.Value;
        
        // Hydrogen is white
        FLinearColor HColor = FLinearColor::White;
        
        auto* HSphere = NewObject<UStaticMeshComponent>(this);
        HSphere->SetStaticMesh(SphereMeshAsset);
        HSphere->SetWorldScale3D(FVector(0.3f)); // Smaller than regular atoms
        HSphere->SetCollisionEnabled(ECollisionEnabled::NoCollision);
        
        auto* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, HSphere);
        Mat->SetVectorParameterValue(TEXT("Color"), HColor);
        Mat->SetScalarParameterValue(TEXT("EmissiveIntensity"), 3.0f);
        HSphere->SetMaterial(0, Mat);
        
        HSphere->AttachToComponent(GetRootComponent(), FAttachmentTransformRules::KeepWorldTransform);
        HSphere->SetWorldLocation(HPos);
        HSphere->RegisterComponent();
        
        HydrogenMeshes.Add(HSphere);
        
        // Draw bond between hydrogen and parent atom
        if (LigInfo->AtomPositions.IsValidIndex(ParentIdx))
        {
            FVector ParentPos = LigInfo->AtomPositions[ParentIdx];
            FString ParentElement = LigInfo->AtomElements.IsValidIndex(ParentIdx) 
                ? LigInfo->AtomElements[ParentIdx] 
                : TEXT("C");
            
            DrawBond(ParentPos, HPos, 1, ParentElement, TEXT("H"), 
                    GetRootComponent(), HydrogenMeshes);
        }
    }
    
    return Hydrogens.Num();
}


void APDBViewer::ClearHydrogens()
{
    for (auto* Mesh : HydrogenMeshes)
    {
        if (Mesh && IsValid(Mesh))
            Mesh->DestroyComponent();
    }
    HydrogenMeshes.Empty();
}



void APDBViewer::DebugPrintLigandInfo()
{
    UE_LOG(LogTemp, Log, TEXT("=== Ligand Debug Info ==="));
    UE_LOG(LogTemp, Log, TEXT("Total ligands: %d"), LigandMap.Num());
    
    for (auto& Pair : LigandMap)
    {
        FLigandInfo* Info = Pair.Value;
        if (!Info)
        {
            UE_LOG(LogTemp, Warning, TEXT("Ligand '%s': NULL INFO"), *Pair.Key);
            continue;
        }
        
        UE_LOG(LogTemp, Log, TEXT("Ligand '%s':"), *Pair.Key);
        UE_LOG(LogTemp, Log, TEXT("  Name: %s"), *Info->LigandName);
        UE_LOG(LogTemp, Log, TEXT("  Visible: %s"), Info->bIsVisible ? TEXT("YES") : TEXT("NO"));
        UE_LOG(LogTemp, Log, TEXT("  Atom Positions: %d"), Info->AtomPositions.Num());
        UE_LOG(LogTemp, Log, TEXT("  Atom Elements: %d"), Info->AtomElements.Num());
        UE_LOG(LogTemp, Log, TEXT("  Bond Pairs: %d"), Info->BondPairs.Num());
        UE_LOG(LogTemp, Log, TEXT("  Bond Orders: %d"), Info->BondOrders.Num());
        UE_LOG(LogTemp, Log, TEXT("  Atom Meshes: %d"), Info->AtomMeshes.Num());
        UE_LOG(LogTemp, Log, TEXT("  Bond Meshes: %d"), Info->BondMeshes.Num());
        
        // Print element breakdown
        TMap<FString, int32> ElementCounts;
        for (const FString& Elem : Info->AtomElements)
        {
            int32& Count = ElementCounts.FindOrAdd(Elem, 0);
            Count++;
        }
        
        UE_LOG(LogTemp, Log, TEXT("  Elements:"));
        for (auto& ElemPair : ElementCounts)
        {
            UE_LOG(LogTemp, Log, TEXT("    %s: %d"), *ElemPair.Key, ElemPair.Value);
        }
        
        // Print bond order breakdown
        if (Info->BondOrders.Num() > 0)
        {
            int32 Single = 0, Double = 0, Triple = 0, Aromatic = 0;
            for (int32 Order : Info->BondOrders)
            {
                if (Order == 1) Single++;
                else if (Order == 2) Double++;
                else if (Order == 3) Triple++;
                else Aromatic++;
            }
            UE_LOG(LogTemp, Log, TEXT("  Bond Types: Single=%d, Double=%d, Triple=%d, Other=%d"), 
                Single, Double, Triple, Aromatic);
        }
    }
    
    UE_LOG(LogTemp, Log, TEXT("Hydrogen meshes: %d"), HydrogenMeshes.Num());
    UE_LOG(LogTemp, Log, TEXT("Hydrogens visible: %s"), bHydrogensVisible ? TEXT("YES") : TEXT("NO"));
    UE_LOG(LogTemp, Log, TEXT("======================"));
}



int32 APDBViewer::GetHydrogenCount() const
{
    return HydrogenMeshes.Num();
}