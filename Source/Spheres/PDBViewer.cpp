// Optimized PDBViewer.cpp

#include "PDBViewer.h"
#include "PDBCameraComponent.h"
#include "HttpModule.h"
#include "Interfaces/IHttpResponse.h"
#include "DrawDebugHelpers.h"
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
#include "Async/Async.h"

// Constants
namespace PDBConstants
{
    constexpr float ATOM_SCALE = 50.0f;
    constexpr float SPHERE_SCALE = 0.5f;
    constexpr float CYLINDER_SCALE = 0.1f;
    constexpr float BOND_OFFSET = 8.0f;
    constexpr float CYLINDER_HALF_HEIGHT = 50.0f;
    constexpr int32 RESERVE_ATOMS = 1000;
    constexpr int32 RESERVE_BONDS = 1000;
}

APDBViewer::APDBViewer()
{
    PrimaryActorTick.bCanEverTick = false;

    USceneComponent* RootComp = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
    RootComponent = RootComp;

    // Load assets with error checking
    static ConstructorHelpers::FObjectFinder<UStaticMesh> SphereMesh(TEXT("/Engine/BasicShapes/Sphere"));
    SphereMeshAsset = SphereMesh.Object;

    static ConstructorHelpers::FObjectFinder<UStaticMesh> CylinderMesh(TEXT("/Engine/BasicShapes/Cylinder"));
    CylinderMeshAsset = CylinderMesh.Object;

    static ConstructorHelpers::FObjectFinder<UMaterial> SphereMat(TEXT("/Engine/BasicShapes/BasicShapeMaterial"));
    SphereMaterialAsset = SphereMat.Object;

    // Reserve memory to reduce allocations
    AllAtomMeshes.Reserve(PDBConstants::RESERVE_ATOMS);
    AllBondMeshes.Reserve(PDBConstants::RESERVE_BONDS);
}

void APDBViewer::BeginPlay()
{
    Super::BeginPlay();
    
    FetchAndDisplayStructure(TEXT("5ENB"));
    
    if (APDBCameraComponent* CameraActor = GetWorld()->SpawnActor<APDBCameraComponent>(
        APDBCameraComponent::StaticClass(),
        GetActorLocation(),
        FRotator::ZeroRotator))
    {
        CameraActor->SetTargetActor(this);
    }
}

// ----------------------------
// Fetch Structure (Optimized)
// ----------------------------
void APDBViewer::FetchAndDisplayStructure(const FString& PDB_ID)
{
    CurrentStructureID = PDB_ID;
    const FString PDB_URL = FString::Printf(TEXT("https://files.rcsb.org/download/%s.pdb"), *PDB_ID);
    const FString CIF_URL = FString::Printf(TEXT("https://files.rcsb.org/download/%s.cif"), *PDB_ID);

    UE_LOG(LogTemp, Log, TEXT("Fetching PDB: %s"), *PDB_URL);

    FetchFileAsync(PDB_URL, [this, CIF_URL](bool bSuccess, const FString& Content)
    {
        if (bSuccess)
        {
            ParsePDB(Content);
        }
        else
        {
            UE_LOG(LogTemp, Warning, TEXT("PDB not found, trying mmCIF"));
            FetchFileAsync(CIF_URL, [this](bool bSuccess2, const FString& Content2)
            {
                if (bSuccess2)
                {
                    ParseMMCIF(Content2);
                }
                else
                {
                    UE_LOG(LogTemp, Error, TEXT("Failed to fetch both PDB and mmCIF"));
                }
            });
        }
    });
}

// ----------------------------
// Async HTTP (Unchanged)
// ----------------------------
void APDBViewer::FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> Callback)
{
    TSharedRef<IHttpRequest, ESPMode::ThreadSafe> Request = FHttpModule::Get().CreateRequest();
    Request->SetURL(URL);
    Request->SetVerb(TEXT("GET"));

    Request->OnProcessRequestComplete().BindLambda(
        [Callback](FHttpRequestPtr Req, FHttpResponsePtr Response, bool bWasSuccessful)
        {
            if (bWasSuccessful && Response.IsValid() && Response->GetResponseCode() == 200)
            {
                Callback(true, Response->GetContentAsString());
            }
            else
            {
                Callback(false, TEXT(""));
            }
        });

    Request->ProcessRequest();
}

// ----------------------------
// Parse PDB (Optimized)
// ----------------------------
void APDBViewer::ParsePDB(const FString& FileContent)
{
    CurrentPDBContent = FileContent;
    ClearResidueMap();

    TArray<FString> Lines;
    FileContent.ParseIntoArrayLines(Lines);

    // Pre-allocate with estimated size
    TMap<FString, TMap<FString, FVector>> ResidueAtoms;
    ResidueAtoms.Reserve(Lines.Num() / 20); // Rough estimate

    TMap<FString, FResidueMetadata> ResidueMetadata;
    ResidueMetadata.Reserve(Lines.Num() / 20);

    int32 AtomCount = 0;
    FString FirstChain;
    bool bHaveFirstChain = false;

    // Single pass parsing
    for (const FString& Line : Lines)
    {
        if (Line.Len() < 80 || !Line.StartsWith(TEXT("HETATM")))
            continue;

        const FString Chain = Line.Mid(21, 1);

        if (!bHaveFirstChain)
        {
            FirstChain = Chain;
            bHaveFirstChain = true;
        }

        if (Chain != FirstChain)
            continue;

        // Parse atom data
        const FString AtomName = Line.Mid(12, 4).TrimStartAndEnd();
        const FString ResidueName = Line.Mid(17, 3).TrimStartAndEnd();
        const FString ResSeq = Line.Mid(22, 4).TrimStartAndEnd();
        
        FString Element = Line.Mid(76, 2).TrimStartAndEnd().ToUpper();
        if (Element.IsEmpty() && !AtomName.IsEmpty())
        {
            Element = AtomName.Left(1).ToUpper();
        }

        const float X = FCString::Atof(*Line.Mid(30, 8)) * PDBConstants::ATOM_SCALE;
        const float Y = FCString::Atof(*Line.Mid(38, 8)) * PDBConstants::ATOM_SCALE;
        const float Z = FCString::Atof(*Line.Mid(46, 8)) * PDBConstants::ATOM_SCALE;

        const FString ResidueKey = FString::Printf(TEXT("%s_%s_%s"), *ResidueName, *ResSeq, *Chain);
        
        ResidueAtoms.FindOrAdd(ResidueKey).Add(AtomName, FVector(X, Y, Z));
        
        // Store metadata only once per residue
        if (!ResidueMetadata.Contains(ResidueKey))
        {
            FResidueMetadata& Meta = ResidueMetadata.Add(ResidueKey);
            Meta.ResidueName = ResidueName;
            Meta.ResidueSeq = ResSeq;
            Meta.Chain = Chain;
            Meta.RecordType = TEXT("HETATM");
        }

        ++AtomCount;
    }

    if (AtomCount == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("No atoms found in PDB"));
        return;
    }

    // Batch create residues
    CreateResiduesFromAtomData(ResidueAtoms, ResidueMetadata);

    UE_LOG(LogTemp, Log, TEXT("Parsed %d atoms from PDB"), AtomCount);
    OnResiduesLoaded.Broadcast();
}

// ----------------------------
// Parse mmCIF (Optimized)
// ----------------------------
void APDBViewer::ParseMMCIF(const FString& FileContent)
{
    CurrentPDBContent = FileContent;
    ClearResidueMap();

    TArray<FString> Lines;
    FileContent.ParseIntoArrayLines(Lines);

    // Parse headers
    TArray<FString> Headers;
    TArray<TArray<FString>> AtomTable;
    AtomTable.Reserve(Lines.Num() / 2);

    int32 XIdx = -1, YIdx = -1, ZIdx = -1, ResIdx = -1, AtomIdx = -1;
    int32 ElementIdx = -1, GroupIdx = -1, ChainIdx = -1, SeqIdx = -1;
    bool bInLoop = false;

    for (const FString& Line : Lines)
    {
        if (Line.StartsWith(TEXT("loop_")))
        {
            bInLoop = true;
            Headers.Empty();
            continue;
        }

        if (bInLoop && Line.StartsWith(TEXT("_atom_site.")))
        {
            const int32 Idx = Headers.Add(Line);
            
            if (Line.Contains(TEXT("Cartn_x"))) XIdx = Idx;
            else if (Line.Contains(TEXT("Cartn_y"))) YIdx = Idx;
            else if (Line.Contains(TEXT("Cartn_z"))) ZIdx = Idx;
            else if (Line.Contains(TEXT("label_comp_id"))) ResIdx = Idx;
            else if (Line.Contains(TEXT("label_atom_id"))) AtomIdx = Idx;
            else if (Line.Contains(TEXT("type_symbol"))) ElementIdx = Idx;
            else if (Line.Contains(TEXT("group_PDB"))) GroupIdx = Idx;
            else if (Line.Contains(TEXT("label_asym_id"))) ChainIdx = Idx;
            else if (Line.Contains(TEXT("label_seq_id")) || Line.Contains(TEXT("auth_seq_id"))) SeqIdx = Idx;
            
            continue;
        }

        if (bInLoop && !Line.StartsWith(TEXT("_")))
        {
            TArray<FString> Tokens;
            Line.ParseIntoArrayWS(Tokens);
            
            if (Tokens.Num() > FMath::Max3(XIdx, YIdx, ZIdx))
            {
                AtomTable.Add(MoveTemp(Tokens));
            }
        }
    }

    // Process atoms
    TMap<FString, TMap<FString, FVector>> ResidueAtoms;
    TMap<FString, FResidueMetadata> ResidueMetadata;
    ResidueAtoms.Reserve(AtomTable.Num() / 20);

    FString FirstChain;
    bool bHaveFirstChain = false;
    int32 AtomCount = 0;

    for (const TArray<FString>& Row : AtomTable)
    {
        if (XIdx < 0 || YIdx < 0 || ZIdx < 0 || ResIdx < 0 || AtomIdx < 0)
            continue;

        if (GroupIdx >= 0 && Row.IsValidIndex(GroupIdx) && Row[GroupIdx] != TEXT("ATOM"))
            continue;

        const FString Chain = (ChainIdx >= 0 && Row.IsValidIndex(ChainIdx)) ? Row[ChainIdx] : TEXT("");
        
        if (!bHaveFirstChain)
        {
            FirstChain = Chain;
            bHaveFirstChain = true;
        }
        
        if (Chain != FirstChain)
            continue;

        const float X = FCString::Atof(*Row[XIdx]) * PDBConstants::ATOM_SCALE;
        const float Y = FCString::Atof(*Row[YIdx]) * PDBConstants::ATOM_SCALE;
        const float Z = FCString::Atof(*Row[ZIdx]) * PDBConstants::ATOM_SCALE;
        
        const FString Residue = Row[ResIdx];
        const FString AtomName = Row[AtomIdx];
        const FString Seq = (SeqIdx >= 0 && Row.IsValidIndex(SeqIdx)) ? Row[SeqIdx] : TEXT("0");
        
        const FString ResidueKey = FString::Printf(TEXT("%s_%s_%s"), *Residue, *Seq, *Chain);
        
        ResidueAtoms.FindOrAdd(ResidueKey).Add(AtomName, FVector(X, Y, Z));
        
        if (!ResidueMetadata.Contains(ResidueKey))
        {
            FResidueMetadata& Meta = ResidueMetadata.Add(ResidueKey);
            Meta.ResidueName = Residue;
            Meta.ResidueSeq = Seq;
            Meta.Chain = Chain;
            Meta.RecordType = TEXT("ATOM");
        }

        ++AtomCount;
    }

    if (AtomCount == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("No atoms found in mmCIF"));
        return;
    }

    CreateResiduesFromAtomData(ResidueAtoms, ResidueMetadata);

    UE_LOG(LogTemp, Log, TEXT("Parsed %d atoms from mmCIF"), AtomCount);
    OnResiduesLoaded.Broadcast();
}

// ----------------------------
// Batch Residue Creation (New)
// ----------------------------
void APDBViewer::CreateResiduesFromAtomData(
    const TMap<FString, TMap<FString, FVector>>& ResidueAtoms,
    const TMap<FString, FResidueMetadata>& Metadata)
{
    for (const auto& Pair : ResidueAtoms)
    {
        const FString& ResidueKey = Pair.Key;
        const TMap<FString, FVector>& Atoms = Pair.Value;
        
        const FResidueMetadata* Meta = Metadata.Find(ResidueKey);
        if (!Meta)
            continue;

        FResidueInfo* ResInfo = new FResidueInfo();
        ResInfo->ResidueName = Meta->ResidueName;
        ResInfo->ResidueSeq = Meta->ResidueSeq;
        ResInfo->Chain = Meta->Chain;
        ResInfo->RecordType = Meta->RecordType;
        ResInfo->bIsVisible = true;
        
        ResInfo->AtomMeshes.Reserve(Atoms.Num());

        // Draw atoms
        for (const auto& Atom : Atoms)
        {
            DrawSphere(Atom.Value.X, Atom.Value.Y, Atom.Value.Z,
                      ElementColor(Atom.Key.Left(1)),
                      GetRootComponent(),
                      ResInfo->AtomMeshes);
        }

        ResidueMap.Add(ResidueKey, ResInfo);

        // Fetch bonds asynchronously
        FetchLigandBondsForResidue(ResidueKey, Meta->ResidueName, Atoms);
    }
}

// ----------------------------
// Fetch Bonds (Optimized)
// ----------------------------
void APDBViewer::FetchLigandBondsForResidue(
    const FString& ResidueKey,
    const FString& ResidueName,
    const TMap<FString, FVector>& AtomPositions)
{
    const FString URL = FString::Printf(
        TEXT("https://files.rcsb.org/ligands/download/%s.cif"),
        *ResidueName.ToUpper());

    TMap<FString, FVector> AtomPosCopy = AtomPositions;
    
    FetchFileAsync(URL, [this, ResidueKey, AtomPosCopy](bool bSuccess, const FString& Content)
    {
        if (bSuccess)
        {
            FResidueInfo** ResInfoPtr = ResidueMap.Find(ResidueKey);
            if (ResInfoPtr && *ResInfoPtr)
            {
                ParseLigandCIFForResidue(Content, AtomPosCopy, *ResInfoPtr);
            }
        }
    });
}

// ----------------------------
// Parse CIF Bonds (Optimized)
// ----------------------------
void APDBViewer::ParseLigandCIFForResidue(
    const FString& FileContent,
    const TMap<FString, FVector>& AtomPositions,
    FResidueInfo* ResInfo)
{
    if (!ResInfo)
        return;

    // Normalize atom IDs once
    TMap<FString, FVector> NormPositions;
    NormPositions.Reserve(AtomPositions.Num());
    
    for (const auto& Pair : AtomPositions)
    {
        const FString Key = NormalizeAtomID(Pair.Key);
        if (!Key.IsEmpty())
        {
            NormPositions.Add(Key, Pair.Value);
        }
    }

    // Parse bonds
    TArray<FString> Lines;
    FileContent.ParseIntoArrayLines(Lines);

    bool bInBondLoop = false;
    TArray<FString> Headers;
    int32 Atom1Idx = -1, Atom2Idx = -1, BondOrderIdx = -1;

    for (const FString& Line : Lines)
    {
        if (Line.StartsWith(TEXT("loop_")))
        {
            bInBondLoop = true;
            Headers.Empty();
            Atom1Idx = Atom2Idx = BondOrderIdx = -1;
            continue;
        }

        if (bInBondLoop && Line.StartsWith(TEXT("_")))
        {
            const int32 Idx = Headers.Add(Line);
            const FString Lower = Line.ToLower();
            
            if (Lower.Contains(TEXT("atom_id_1")) || Lower.Contains(TEXT("atom_1")))
                Atom1Idx = Idx;
            else if (Lower.Contains(TEXT("atom_id_2")) || Lower.Contains(TEXT("atom_2")))
                Atom2Idx = Idx;
            else if (Lower.Contains(TEXT("value_order")) || Lower.Contains(TEXT("bond_order")))
                BondOrderIdx = Idx;
            
            continue;
        }

        if (bInBondLoop && !Line.StartsWith(TEXT("_")))
        {
            if (Atom1Idx < 0 || Atom2Idx < 0)
                continue;

            TArray<FString> Tokens;
            Line.ParseIntoArrayWS(Tokens);

            if (Tokens.Num() <= FMath::Max(Atom1Idx, Atom2Idx))
                continue;

            const FString Id1 = NormalizeAtomID(Tokens[Atom1Idx]);
            const FString Id2 = NormalizeAtomID(Tokens[Atom2Idx]);
            const int32 Order = ParseBondOrder(
                BondOrderIdx >= 0 && Tokens.IsValidIndex(BondOrderIdx) ? Tokens[BondOrderIdx] : TEXT("1"));

            const FVector* P1 = NormPositions.Find(Id1);
            const FVector* P2 = NormPositions.Find(Id2);

            if (P1 && P2)
            {
                DrawBond(*P1, *P2, Order, FLinearColor::Gray, GetRootComponent(), ResInfo->BondMeshes);
            }
        }

        if (bInBondLoop && Line.StartsWith(TEXT("data_")))
        {
            break; // Exit loop
        }
    }
}

// ----------------------------
// Helper: Normalize Atom ID
// ----------------------------
FString APDBViewer::NormalizeAtomID(const FString& In) const
{
    FString Result;
    Result.Reserve(In.Len());
    
    for (const TCHAR& C : In)
    {
        if (FChar::IsAlnum(C))
        {
            Result.AppendChar(FChar::ToUpper(C));
        }
    }
    
    return Result;
}

// ----------------------------
// Helper: Parse Bond Order
// ----------------------------
int32 APDBViewer::ParseBondOrder(const FString& OrderStr) const
{
    const FString Lower = OrderStr.ToLower();
    
    if (Lower.Contains(TEXT("ar")) || Lower.Contains(TEXT("aromatic")))
        return 1;
    
    if (Lower.Contains(TEXT("doub")) || Lower.Equals(TEXT("d")))
        return 2;
    
    if (Lower.Contains(TEXT("trip")) || Lower.Equals(TEXT("t")))
        return 3;
    
    const int32 Parsed = FCString::Atoi(*OrderStr);
    return (Parsed > 0 && Parsed <= 3) ? Parsed : 1;
}

// ----------------------------
// Draw Sphere (Optimized)
// ----------------------------
void APDBViewer::DrawSphere(
    float X, float Y, float Z,
    const FLinearColor& Color,
    USceneComponent* Parent,
    TArray<UStaticMeshComponent*>& OutArray)
{
    if (!SphereMeshAsset || !SphereMaterialAsset || !Parent)
        return;

    UStaticMeshComponent* Sphere = NewObject<UStaticMeshComponent>(this);
    Sphere->SetStaticMesh(SphereMeshAsset);
    Sphere->SetWorldScale3D(FVector(PDBConstants::SPHERE_SCALE));
    Sphere->SetCollisionEnabled(ECollisionEnabled::NoCollision);
    
    // Create material instance
    UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
    Mat->SetVectorParameterValue(TEXT("Color"), Color);
    Mat->SetScalarParameterValue(TEXT("EmissiveIntensity"), 5.0f);
    Sphere->SetMaterial(0, Mat);

    Sphere->AttachToComponent(Parent, FAttachmentTransformRules::KeepWorldTransform);
    
    const FVector WorldPos(X, Y, Z);
    const FVector LocalPos = Parent->GetComponentTransform().InverseTransformPosition(WorldPos);
    Sphere->SetRelativeLocation(LocalPos);
    
    Sphere->RegisterComponent();

    OutArray.Add(Sphere);
    AllAtomMeshes.Add(Sphere);
}

// ----------------------------
// Draw Bond (Optimized)
// ----------------------------
void APDBViewer::DrawBond(
    const FVector& Start,
    const FVector& End,
    int32 Order,
    const FLinearColor& Color,
    USceneComponent* Parent,
    TArray<UStaticMeshComponent*>& OutArray)
{
    if (!CylinderMeshAsset || !SphereMaterialAsset || !Parent)
        return;

    const FVector BondVector = End - Start;
    const float Length = BondVector.Size();
    const FVector MidPoint = Start + 0.5f * BondVector;
    const FRotator Rotation = FRotationMatrix::MakeFromZ(BondVector).Rotator();
    const float ZScale = Length / (2.0f * PDBConstants::CYLINDER_HALF_HEIGHT);

    auto CreateCylinder = [&](const FVector& Position)
    {
        UStaticMeshComponent* Cylinder = NewObject<UStaticMeshComponent>(this);
        Cylinder->SetStaticMesh(CylinderMeshAsset);
        Cylinder->SetWorldLocation(Position);
        Cylinder->SetWorldRotation(Rotation);
        Cylinder->SetWorldScale3D(FVector(PDBConstants::CYLINDER_SCALE, PDBConstants::CYLINDER_SCALE, ZScale));
        Cylinder->SetCollisionEnabled(ECollisionEnabled::NoCollision);
        
        UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
        Mat->SetVectorParameterValue(TEXT("Color"), Color);
        Cylinder->SetMaterial(0, Mat);
        
        Cylinder->AttachToComponent(Parent, FAttachmentTransformRules::KeepWorldTransform);
        Cylinder->RegisterComponent();
        
        OutArray.Add(Cylinder);
        AllBondMeshes.Add(Cylinder);
    };

    if (Order <= 1)
    {
        CreateCylinder(MidPoint);
        return;
    }

    // Calculate perpendicular vector
    FVector Perp = FVector::CrossProduct(BondVector.GetSafeNormal(), FVector::UpVector);
    if (Perp.SizeSquared() < KINDA_SMALL_NUMBER)
    {
        Perp = FVector::CrossProduct(BondVector.GetSafeNormal(), FVector::RightVector);
    }
    Perp.Normalize();

    if (Order == 2)
    {
        CreateCylinder(MidPoint + Perp * PDBConstants::BOND_OFFSET);
        CreateCylinder(MidPoint - Perp * PDBConstants::BOND_OFFSET);
    }
    else if (Order == 3)
    {
        CreateCylinder(MidPoint);
        CreateCylinder(MidPoint + Perp * PDBConstants::BOND_OFFSET);
        CreateCylinder(MidPoint - Perp * PDBConstants::BOND_OFFSET);
    }
    else
    {
        CreateCylinder(MidPoint);
    }
}

// ----------------------------
// Element Colors (Optimized with static map)
// ----------------------------
FLinearColor APDBViewer::ElementColor(const FString& Element) const
{
    static const TMap<FString, FLinearColor> ColorMap = {
        {TEXT("C"), FLinearColor(0.1f, 0.1f, 0.1f)},
        {TEXT("O"), FLinearColor::Red},
        {TEXT("H"), FLinearColor::White},
        {TEXT("D"), FLinearColor::White},
        {TEXT("N"), FLinearColor::Blue},
        {TEXT("S"), FLinearColor::Yellow},
        {TEXT("CL"), FLinearColor(0.0f, 1.0f, 0.0f)},
        {TEXT("P"), FLinearColor(1.0f, 0.5f, 0.0f)},
        {TEXT("F"), FLinearColor(0.0f, 1.0f, 0.0f)},
        {TEXT("BR"), FLinearColor(0.6f, 0.2f, 0.2f)},
        {TEXT("I"), FLinearColor(0.4f, 0.0f, 0.8f)},
        {TEXT("FE"), FLinearColor(0.8f, 0.4f, 0.0f)},
        {TEXT("MG"), FLinearColor(0.0f, 0.8f, 0.0f)},
        {TEXT("ZN"), FLinearColor(0.5f, 0.5f, 0.5f)},
        {TEXT("CA"), FLinearColor(0.2f, 0.6f, 1.0f)},
        {TEXT("NA"), FLinearColor(0.0f, 0.0f, 1.0f)},
        {TEXT("K"), FLinearColor(0.5f, 0.0f, 1.0f)},
        {TEXT("CU"), FLinearColor(1.0f, 0.5f, 0.0f)},
        {TEXT("B"), FLinearColor(1.0f, 0.7f, 0.7f)}
    };

    const FLinearColor* Found = ColorMap.Find(Element.ToUpper());
    return Found ? *Found : FLinearColor::Gray;
}

// ----------------------------
// Cleanup
// ----------------------------
void APDBViewer::ClearResidueMap()
{
    for (auto& Pair : ResidueMap)
    {
        delete Pair.Value;
    }
    ResidueMap.Empty();
}

void APDBViewer::ClearCurrentStructure()
{
    ClearResidueMap();

    for (UStaticMeshComponent* Mesh : AllAtomMeshes)
    {
        if (Mesh && IsValid(Mesh))
        {
            Mesh->DestroyComponent();
        }
    }
    AllAtomMeshes.Empty();

    for (UStaticMeshComponent* Mesh : AllBondMeshes)
    {
        if (Mesh && IsValid(Mesh))
        {
            Mesh->DestroyComponent();
        }
    }
    AllBondMeshes.Empty();

    UE_LOG(LogTemp, Log, TEXT("Cleared current structure"));
}

// ----------------------------
// Residue Management
// ----------------------------
void APDBViewer::ToggleResidueVisibility(const FString& ResidueKey)
{
    FResidueInfo** ResInfoPtr = ResidueMap.Find(ResidueKey);
    if (!ResInfoPtr || !*ResInfoPtr)
        return;

    FResidueInfo* ResInfo = *ResInfoPtr;
    ResInfo->bIsVisible = !ResInfo->bIsVisible;

    for (UStaticMeshComponent* Mesh : ResInfo->AtomMeshes)
    {
        if (Mesh && IsValid(Mesh))
        {
            Mesh->SetVisibility(ResInfo->bIsVisible);
        }
    }

    for (UStaticMeshComponent* Mesh : ResInfo->BondMeshes)
    {
        if (Mesh && IsValid(Mesh))
        {
            Mesh->SetVisibility(ResInfo->bIsVisible);
        }
    }

    UE_LOG(LogTemp, Log, TEXT("Toggled %s: %s"), *ResidueKey, ResInfo->bIsVisible ? TEXT("ON") : TEXT("OFF"));
}

TArray<FString> APDBViewer::GetResidueList() const
{
    TArray<FString> List;
    ResidueMap.GetKeys(List);
    
    List.Sort([](const FString& A, const FString& B)
    {
        int32 UnderA = A.Find(TEXT("_"));
        int32 UnderB = B.Find(TEXT("_"));
        
        if (UnderA != INDEX_NONE && UnderB != INDEX_NONE)
        {
            FString SeqA = A.Mid(UnderA + 1);
            FString SeqB = B.Mid(UnderB + 1);
            
            int32 SecondUnderA = SeqA.Find(TEXT("_"));
            int32 SecondUnderB = SeqB.Find(TEXT("_"));
            
            if (SecondUnderA != INDEX_NONE) SeqA = SeqA.Left(SecondUnderA);
            if (SecondUnderB != INDEX_NONE) SeqB = SeqB.Left(SecondUnderB);
            
            return FCString::Atoi(*SeqA) < FCString::Atoi(*SeqB);
        }
        return A < B;
    });

    return List;
}

TArray<FString> APDBViewer::GetLigandList() const
{
    return GetResidueList(); // Same implementation
}

FString APDBViewer::GetResidueDisplayName(const FString& ResidueKey) const
{
    FResidueInfo* const* ResInfoPtr = ResidueMap.Find(ResidueKey);
    if (!ResInfoPtr || !*ResInfoPtr)
        return ResidueKey;
    
    FResidueInfo* ResInfo = *ResInfoPtr;
    return FString::Printf(TEXT("%s %s"), *ResInfo->ResidueName, *ResInfo->ResidueSeq);
}

// ----------------------------
// File Operations
// ----------------------------
void APDBViewer::SaveStructureToFile(const FString& FilePath)
{
    if (CurrentPDBContent.IsEmpty())
    {
        UE_LOG(LogTemp, Warning, TEXT("No structure data to save"));
        return;
    }

    if (FFileHelper::SaveStringToFile(CurrentPDBContent, *FilePath))
    {
        UE_LOG(LogTemp, Log, TEXT("Saved structure to: %s"), *FilePath);
    }
    else
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to save structure to: %s"), *FilePath);
    }
}

void APDBViewer::LoadStructureFromFile(const FString& FilePath)
{
    FString FileContent;
    if (!FFileHelper::LoadFileToString(FileContent, *FilePath))
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to load file: %s"), *FilePath);
        return;
    }

    UE_LOG(LogTemp, Log, TEXT("Loading structure from: %s"), *FilePath);
    ClearCurrentStructure();
    
    const FString Extension = FPaths::GetExtension(FilePath).ToLower();
    CurrentStructureID = FPaths::GetBaseFilename(FilePath);
    CurrentPDBContent = FileContent;

    if (Extension == TEXT("pdb"))
    {
        ParsePDB(FileContent);
    }
    else if (Extension == TEXT("cif"))
    {
        ParseMMCIF(FileContent);
    }
    else
    {
        UE_LOG(LogTemp, Error, TEXT("Unsupported file format: %s"), *Extension);
    }
}

bool APDBViewer::ShowFileDialog(bool bSave, FString& OutFilePath)
{
    IDesktopPlatform* DesktopPlatform = FDesktopPlatformModule::Get();
    if (!DesktopPlatform)
    {
        UE_LOG(LogTemp, Error, TEXT("Desktop Platform not available"));
        return false;
    }

    TSharedPtr<SWindow> ParentWindow;
    if (FModuleManager::Get().IsModuleLoaded("MainFrame"))
    {
        IMainFrameModule& MainFrame = FModuleManager::LoadModuleChecked<IMainFrameModule>("MainFrame");
        ParentWindow = MainFrame.GetParentWindow();
    }

    void* ParentWindowHandle = (ParentWindow.IsValid() && ParentWindow->GetNativeWindow().IsValid())
        ? ParentWindow->GetNativeWindow()->GetOSWindowHandle()
        : nullptr;

    TArray<FString> OutFiles;
    const FString FileTypes = TEXT("PDB Files (*.pdb)|*.pdb|mmCIF Files (*.cif)|*.cif|All Files (*.*)|*.*");
    const FString DefaultPath = FPaths::ProjectSavedDir();

    bool bSuccess = false;
    if (bSave)
    {
        bSuccess = DesktopPlatform->SaveFileDialog(
            ParentWindowHandle,
            TEXT("Save PDB Structure"),
            DefaultPath,
            TEXT("structure.pdb"),
            FileTypes,
            EFileDialogFlags::None,
            OutFiles);
    }
    else
    {
        bSuccess = DesktopPlatform->OpenFileDialog(
            ParentWindowHandle,
            TEXT("Load PDB Structure"),
            DefaultPath,
            TEXT(""),
            FileTypes,
            EFileDialogFlags::None,
            OutFiles);
    }

    if (bSuccess && OutFiles.Num() > 0)
    {
        OutFilePath = OutFiles[0];
        return true;
    }

    return false;
}

void APDBViewer::OpenSaveDialog()
{
    FString FilePath;
    if (ShowFileDialog(true, FilePath))
    {
        SaveStructureToFile(FilePath);
    }
}

void APDBViewer::OpenLoadDialog()
{
    FString FilePath;
    if (ShowFileDialog(false, FilePath))
    {
        LoadStructureFromFile(FilePath);
    }
}