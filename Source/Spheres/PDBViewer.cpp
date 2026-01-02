// Add these includes at the top of your PDBViewer.cpp file:
#include "Misc/FileHelper.h"
#include "HAL/PlatformFilemanager.h"
#include "Misc/Paths.h"
#include "IDesktopPlatform.h"
#include "DesktopPlatformModule.h"
#include "Interfaces/IMainFrameModule.h"

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
#include "UObject/ConstructorHelpers.h"  // ADD THIS

APDBViewer::APDBViewer()
{
    PrimaryActorTick.bCanEverTick = false;

    USceneComponent* RootComp = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
    RootComponent = RootComp;

    static ConstructorHelpers::FObjectFinder<UStaticMesh> SphereMeshObj(TEXT("/Engine/BasicShapes/Sphere.Sphere"));
    if (SphereMeshObj.Succeeded()) SphereMeshAsset = SphereMeshObj.Object;

    static ConstructorHelpers::FObjectFinder<UStaticMesh> CylinderMeshObj(TEXT("/Engine/BasicShapes/Cylinder.Cylinder"));
    if (CylinderMeshObj.Succeeded()) CylinderMeshAsset = CylinderMeshObj.Object;

    // CHANGE: Use UMaterial instead of UMaterialInterface
    static ConstructorHelpers::FObjectFinder<UMaterial> SphereMatObj(TEXT("/Engine/BasicShapes/BasicShapeMaterial.BasicShapeMaterial"));
    if (SphereMatObj.Succeeded()) SphereMaterialAsset = SphereMatObj.Object;
}

void APDBViewer::BeginPlay()
{
    Super::BeginPlay();
    FetchAndDisplayStructure(TEXT("5ENB")); // example PDB ID
    APDBCameraComponent* CameraActor = GetWorld()->SpawnActor<APDBCameraComponent>(
        APDBCameraComponent::StaticClass(), 
        GetActorLocation(), 
        FRotator::ZeroRotator
    );
    if (CameraActor)
    {
        CameraActor->SetTargetActor(this);
    }
}

// ----------------------------
// Fetch PDB or CIF
// ----------------------------
void APDBViewer::FetchAndDisplayStructure(const FString& PDB_ID)
{
    CurrentStructureID = PDB_ID;
    FString PDB_URL = FString::Printf(TEXT("https://files.rcsb.org/download/%s.pdb"), *PDB_ID);
    FString CIF_URL = FString::Printf(TEXT("https://files.rcsb.org/download/%s.cif"), *PDB_ID);

    UE_LOG(LogTemp, Log, TEXT("Fetching PDB: %s"), *PDB_URL);

    FetchFileAsync(PDB_URL, [this, CIF_URL](bool bSuccess, const FString& Content)
    {
        if (bSuccess) ParsePDB(Content);
        else
        {
            UE_LOG(LogTemp, Warning, TEXT("PDB not found, trying mmCIF"));
            FetchFileAsync(CIF_URL, [this](bool bSuccess2, const FString& Content2)
            {
                if (bSuccess2) ParseMMCIF(Content2);
                else UE_LOG(LogTemp, Error, TEXT("Failed to fetch both PDB and mmCIF"));
            });
        }
    });
}

// ----------------------------
// Async HTTP fetch
// ----------------------------
void APDBViewer::FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> Callback)
{
    TSharedRef<IHttpRequest, ESPMode::ThreadSafe> Request = FHttpModule::Get().CreateRequest();
    Request->SetURL(URL);
    Request->SetVerb(TEXT("GET"));

    Request->OnProcessRequestComplete().BindLambda([Callback](FHttpRequestPtr Req, FHttpResponsePtr Response, bool bWasSuccessful)
    {
        if (bWasSuccessful && Response.IsValid() && Response->GetResponseCode() == 200)
            Callback(true, Response->GetContentAsString());
        else
            Callback(false, TEXT(""));
    });

    Request->ProcessRequest();
}

// ----------------------------
// Parse PDB (ATOM only)
// ----------------------------
void APDBViewer::ParsePDB(const FString& FileContent)
{
    // Store the content for saving later
    CurrentPDBContent = FileContent;
    
    TArray<FString> Lines;
    FileContent.ParseIntoArrayLines(Lines);

    TMap<FString, TMap<FString, FVector>> ResidueAtoms;
    TMap<FString, FString> ResidueKeyToName;

    FVector Center(0,0,0);
    int32 AtomCount = 0;
    const float Scale = 50.f;

    bool bHaveFirstChain = false;
    FString FirstChain;

    for (const FString& Line : Lines)
    {
        if (!Line.StartsWith("ATOM")) continue;

        FString Chain = Line.Mid(21, 1);

        if (!bHaveFirstChain)
        {
            FirstChain = Chain;
            bHaveFirstChain = true;
        }

        if (Chain != FirstChain) continue;

        FString AtomName = Line.Mid(12,4).TrimStartAndEnd();
        FString ResidueName = Line.Mid(17,3).TrimStartAndEnd();
        FString ResSeq = Line.Mid(22,4).TrimStartAndEnd();
        FString Element = Line.Mid(76,2).TrimStartAndEnd().ToUpper();
        if (Element.IsEmpty() && AtomName.Len() > 0) Element = AtomName.Left(1).ToUpper();

        float X = FCString::Atof(*Line.Mid(30,8)) * Scale;
        float Y = FCString::Atof(*Line.Mid(38,8)) * Scale;
        float Z = FCString::Atof(*Line.Mid(46,8)) * Scale;

        FVector Pos(X,Y,Z);
        Center += Pos;
        AtomCount++;

        FString ResidueKey = FString::Printf(TEXT("%s_%s_%s"), *ResidueName, *ResSeq, *Chain);
        ResidueAtoms.FindOrAdd(ResidueKey).Add(AtomName, Pos);
        ResidueKeyToName.FindOrAdd(ResidueKey) = ResidueName;
    }

    if (AtomCount == 0) return;

    Center /= AtomCount;

    for (auto& Pair : ResidueAtoms)
    {
        const FString& UniqueResidueKey = Pair.Key;
        TMap<FString, FVector>& Atoms = Pair.Value;

        for (auto& Atom : Atoms)
        {
            FVector Pos = Atom.Value;
            DrawSphere(Pos.X, Pos.Y, Pos.Z, ElementColor(Atom.Key.Left(1)), GetRootComponent(), AllAtomMeshes);
        }

        const FString* Residue3Name = ResidueKeyToName.Find(UniqueResidueKey);
        if (Residue3Name)
            FetchLigandCIF(*Residue3Name, Atoms);
    }

    UE_LOG(LogTemp, Log, TEXT("Parsed %d atoms from PDB (first chain only)"), AtomCount);
}

// ----------------------------
// Parse mmCIF (ATOM only)
// ----------------------------
void APDBViewer::ParseMMCIF(const FString& FileContent)
{
    CurrentPDBContent = FileContent;
    
    TArray<FString> Lines;
    FileContent.ParseIntoArrayLines(Lines);

    TArray<FString> Headers;
    TArray<TArray<FString>> AtomTable;

    int XIdx=-1, YIdx=-1, ZIdx=-1, ResIdx=-1, AtomIdx=-1, ElementIdx=-1, GroupIdx=-1, ChainIdx=-1, SeqIdx=-1;
    bool bInLoop = false;

    for (const FString& Line : Lines)
    {
        if (Line.StartsWith("loop_")) { bInLoop = true; Headers.Empty(); continue; }
        if (bInLoop && Line.StartsWith("_atom_site."))
        {
            Headers.Add(Line);
            if (Line.Contains("Cartn_x")) XIdx = Headers.Num()-1;
            if (Line.Contains("Cartn_y")) YIdx = Headers.Num()-1;
            if (Line.Contains("Cartn_z")) ZIdx = Headers.Num()-1;
            if (Line.Contains("label_comp_id")) ResIdx = Headers.Num()-1;
            if (Line.Contains("label_atom_id")) AtomIdx = Headers.Num()-1;
            if (Line.Contains("type_symbol")) ElementIdx = Headers.Num()-1;
            if (Line.Contains("group_PDB")) GroupIdx = Headers.Num()-1;
            if (Line.Contains("label_asym_id")) ChainIdx = Headers.Num()-1;
            if (Line.Contains("label_seq_id") || Line.Contains("auth_seq_id")) SeqIdx = Headers.Num()-1;
            continue;
        }

        if (bInLoop && !Line.StartsWith("_"))
        {
            TArray<FString> Tokens;
            Line.ParseIntoArrayWS(Tokens);
            if (Tokens.Num() > FMath::Max3(XIdx,YIdx,ZIdx)) AtomTable.Add(Tokens);
        }
    }

    TMap<FString, TMap<FString, FVector>> ResidueAtoms;
    TMap<FString, FString> ResidueKeyToName;
    TArray<UStaticMeshComponent*> AtomMeshes;
    int32 AtomCount=0;

    const float Scale = 50.f; // match BlueSphere scaling

    bool bHaveFirstChain = false;
    FString FirstChain;

    for (const auto& Row : AtomTable)
    {
        if (XIdx<0||YIdx<0||ZIdx<0||ResIdx<0||AtomIdx<0) continue;
        if (GroupIdx>=0 && Row[GroupIdx] != "ATOM") continue;

        // If chain info is present, only accept rows from the first chain encountered
        if (ChainIdx >= 0)
        {
            if (!Row.IsValidIndex(ChainIdx)) continue;
            FString Chain = Row[ChainIdx];
            if (!bHaveFirstChain)
            {
                FirstChain = Chain;
                bHaveFirstChain = true;
            }
            if (Chain != FirstChain) continue;
        }

        float X = FCString::Atof(*Row[XIdx]) * Scale;
        float Y = FCString::Atof(*Row[YIdx]) * Scale;
        float Z = FCString::Atof(*Row[ZIdx]) * Scale;
        FString Residue = Row[ResIdx];
        FString AtomName = Row[AtomIdx];
        FString Element = (ElementIdx>=0 && Row.IsValidIndex(ElementIdx)) ? Row[ElementIdx].ToUpper() : AtomName.Left(1).ToUpper();

        FString Seq = (SeqIdx >= 0 && Row.IsValidIndex(SeqIdx)) ? Row[SeqIdx] : FString(TEXT("0"));
        FString Chain = (ChainIdx >= 0 && Row.IsValidIndex(ChainIdx)) ? Row[ChainIdx] : FString(TEXT(""));

        // Unique residue key that includes seq and chain (prevents collisions)
        FString ResidueKey = FString::Printf(TEXT("%s_%s_%s"), *Residue, *Seq, *Chain);

        FVector Pos(X,Y,Z);
        AtomCount++;

        ResidueAtoms.FindOrAdd(ResidueKey).Add(AtomName, Pos);
        ResidueKeyToName.FindOrAdd(ResidueKey) = Residue;
    }

    if (AtomCount==0) return;

    for (auto& Pair : ResidueAtoms)
    {
        const FString& UniqueResidueKey = Pair.Key;
        TMap<FString, FVector>& Atoms = Pair.Value;

        for (auto& Atom : Atoms)
        {
            FVector Pos = Atom.Value;
            DrawSphere(Pos.X, Pos.Y, Pos.Z, ElementColor(Atom.Key.Left(1)), GetRootComponent(), AllAtomMeshes);
        }

        // Fetch bonds for this residue using 3-letter residue name
        const FString* Residue3Name = ResidueKeyToName.Find(UniqueResidueKey);
        if (Residue3Name)
            FetchLigandCIF(*Residue3Name, Atoms);
    }

    UE_LOG(LogTemp, Log, TEXT("Parsed %d atoms from mmCIF (first chain only)"), AtomCount);
}

// ----------------------------
// Fetch ligand CIF for bonds
// ----------------------------
void APDBViewer::FetchLigandCIF(const FString& ResidueName, const TMap<FString, FVector>& AtomPositions)
{
    FString URL = FString::Printf(TEXT("https://files.rcsb.org/ligands/download/%s.cif"), *ResidueName.ToUpper());
    UE_LOG(LogTemp, Log, TEXT("Fetching ligand CIF: %s"), *URL);

    // Copy AtomPositions to avoid capture issues
    TMap<FString, FVector> AtomPosCopy = AtomPositions;
    FetchFileAsync(URL, [this, AtomPosCopy](bool bSuccess, const FString& Content)
    {
        if (bSuccess) ParseLigandCIF(Content, AtomPosCopy);
    });
}

// ----------------------------
// Parse CIF bonds
// ----------------------------
void APDBViewer::ParseLigandCIF(const FString& FileContent, const TMap<FString, FVector>& AtomPositions)
{
    TArray<FString> Lines;
    FileContent.ParseIntoArrayLines(Lines);

    // Build normalized lookup of provided atom positions (keys => normalized uppercase)
    auto NormalizeId = [](const FString& In) -> FString
    {
        FString S = In;
        S = S.Replace(TEXT("\""), TEXT(""));
        S = S.Replace(TEXT("'"), TEXT(""));
        S.TrimStartAndEndInline();
        // Keep only letters and digits (remove punctuation like commas/slashes/parentheses)
        FString Out;
        for (int32 i = 0; i < S.Len(); ++i)
        {
            TCHAR C = S[i];
            if (FChar::IsAlnum(C))
                Out.AppendChar(C);
        }
        Out = Out.ToUpper();
        return Out;
    };

    TMap<FString, FVector> NormPositions;
    for (const auto& Pair : AtomPositions)
    {
        FString Key = NormalizeId(Pair.Key);
        if (!Key.IsEmpty() && !NormPositions.Contains(Key))
            NormPositions.Add(Key, Pair.Value);
    }

    bool bInBondLoop = false;
    TArray<FString> Headers;
    int Atom1Idx = -1;
    int Atom2Idx = -1;
    int BondOrderIdx = -1; // <-- new

    TArray<UStaticMeshComponent*> BondMeshes;

    for (const FString& Line : Lines)
    {
        if (Line.StartsWith("loop_"))
        {
            bInBondLoop = true;
            Headers.Empty();
            Atom1Idx = Atom2Idx = BondOrderIdx = -1;
            continue;
        }

        if (bInBondLoop && Line.StartsWith("_"))
        {
            // collect header lines for this loop
            Headers.Add(Line);
            FString Lower = Line.ToLower();
            if (Lower.Contains(TEXT("atom_id_1")) || Lower.Contains(TEXT("atom_1")))
                Atom1Idx = Headers.Num() - 1;
            if (Lower.Contains(TEXT("atom_id_2")) || Lower.Contains(TEXT("atom_2")))
                Atom2Idx = Headers.Num() - 1;
            if (Lower.Contains(TEXT("value_order")) || Lower.Contains(TEXT("bond_order")) || Lower.Contains(TEXT("value")))
                BondOrderIdx = Headers.Num() - 1; // <-- new
            continue;
        }

        if (bInBondLoop && !Line.StartsWith("_"))
        {
            // We only parse rows when we already found the required columns
            if (Atom1Idx < 0 || Atom2Idx < 0)
            {
                // if bond loop header didn't contain both ids, skip parsing rows
                continue;
            }

            TArray<FString> Tokens;
            Line.ParseIntoArrayWS(Tokens);

            int32 MaxIdx = FMath::Max(Atom1Idx, Atom2Idx);
            if (Tokens.Num() <= MaxIdx) continue;

            FString Raw1 = Tokens[Atom1Idx];
            FString Raw2 = Tokens[Atom2Idx];
            FString RawOrder = (BondOrderIdx >= 0 && Tokens.IsValidIndex(BondOrderIdx)) ? Tokens[BondOrderIdx] : FString(TEXT("1"));

            int Order = 1;
            FString LowerOrder = RawOrder.ToLower();
            if (LowerOrder.Contains(TEXT("ar")) || LowerOrder.Contains(TEXT("aromatic")))
            {
                Order = 1; // handle aromatic like single for now (could special-case)
            }
            else
            {
                int Parsed = FCString::Atoi(*RawOrder);
                if (Parsed > 0) Order = Parsed;
                else if (LowerOrder.Contains(TEXT("double")) || LowerOrder.Contains(TEXT("d"))) Order = 2;
                else if (LowerOrder.Contains(TEXT("triple")) || LowerOrder.Contains(TEXT("t"))) Order = 3;
            }

            // Normalize the raw atom ids and try to find positions
            FString Id1 = NormalizeId(Raw1);
            FString Id2 = NormalizeId(Raw2);

            const FVector* P1 = NormPositions.Find(Id1);
            const FVector* P2 = NormPositions.Find(Id2);

            // If direct normalized match failed, try a case-insensitive fallback over all keys
            if ((!P1 || !P2) && NormPositions.Num() > 0)
            {
                if (!P1)
                {
                    for (const auto& NP : NormPositions)
                    {
                        if (NP.Key.Equals(Id1, ESearchCase::IgnoreCase))
                        {
                            P1 = &NP.Value;
                            break;
                        }
                    }
                }
                if (!P2)
                {
                    for (const auto& NP : NormPositions)
                    {
                        if (NP.Key.Equals(Id2, ESearchCase::IgnoreCase))
                        {
                            P2 = &NP.Value;
                            break;
                        }
                    }
                }
            }

            if (P1 && P2)
            {
                DrawBond(*P1, *P2, Order, FLinearColor::Gray, GetRootComponent(), AllBondMeshes);
            }
        }

        // Exit loop parsing when we hit a new data block (optional)
        if (bInBondLoop && Line.StartsWith("data_"))
        {
            bInBondLoop = false;
            Headers.Empty();
            Atom1Idx = Atom2Idx = -1;
        }
    }
}

// ----------------------------
// Draw a sphere
// ----------------------------
void APDBViewer::DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray)
{
    if (!SphereMeshAsset || !SphereMaterialAsset || !Parent) return;

    const FVector WorldPos(x, y, z);
    // Convert world-style position into Parent-local position so both spheres and bonds share the same space
    const FTransform ParentTransform = Parent->GetComponentTransform();
    const FVector LocalPos = ParentTransform.InverseTransformPosition(WorldPos);

    UStaticMeshComponent* Sphere = NewObject<UStaticMeshComponent>(this);
    Sphere->SetStaticMesh(SphereMeshAsset);
    Sphere->AttachToComponent(Parent, FAttachmentTransformRules::KeepRelativeTransform);
    Sphere->SetRelativeLocation(LocalPos);
    Sphere->SetWorldScale3D(FVector(0.5f)); // fine to use world scale for uniform sizing

    UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
    Mat->SetVectorParameterValue(FName("Color"), Color);
    Mat->SetScalarParameterValue(FName("EmissiveIntensity"), 5.f);
    Sphere->SetMaterial(0, Mat);

    Sphere->RegisterComponent();

    OutArray.Add(Sphere);
}

// ----------------------------
// Draw a bond
// ----------------------------
void APDBViewer::DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray)
{
    if (!CylinderMeshAsset || !SphereMaterialAsset || !Parent) return;

    FVector BondVector = End - Start;
    float Length = BondVector.Size();
    FVector MidPoint = Start + 0.5f * BondVector;
    FRotator Rotation = FRotationMatrix::MakeFromZ(BondVector).Rotator();

    const float DefaultHalfHeight = 50.0f;
    float ZScale = Length / (2.0f * DefaultHalfHeight);

    auto SpawnCylinderAt = [&](const FVector& Pos){
        UStaticMeshComponent* Cylinder = NewObject<UStaticMeshComponent>(this);
        Cylinder->SetStaticMesh(CylinderMeshAsset);
        Cylinder->SetWorldLocation(Pos);
        Cylinder->SetWorldRotation(Rotation);
        Cylinder->SetWorldScale3D(FVector(0.1f, 0.1f, ZScale));
        Cylinder->SetCollisionEnabled(ECollisionEnabled::NoCollision);
        Cylinder->RegisterComponent();
        Cylinder->AttachToComponent(Parent, FAttachmentTransformRules::KeepWorldTransform);

        UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
        Mat->SetVectorParameterValue(FName("Color"), Color);
        Cylinder->SetMaterial(0, Mat);
        OutArray.Add(Cylinder);
    };

    // If single bond, draw single cylinder
    if (Order <= 1)
    {
        SpawnCylinderAt(MidPoint);
        return;
    }

    // For double/triple bonds: compute a perpendicular offset vector and spawn multiple cylinders
    FVector Dir = BondVector.GetSafeNormal();
    FVector Perp = FVector::CrossProduct(Dir, FVector::UpVector);
    if (Perp.SizeSquared() < KINDA_SMALL_NUMBER)
        Perp = FVector::CrossProduct(Dir, FVector::RightVector);
    Perp.Normalize();

    const float OffsetDist = 8.0f; // tweak for visual spacing

    if (Order == 2)
    {
        SpawnCylinderAt(MidPoint + Perp * OffsetDist);
        SpawnCylinderAt(MidPoint - Perp * OffsetDist);
    }
    else if (Order == 3)
    {
        SpawnCylinderAt(MidPoint);
        SpawnCylinderAt(MidPoint + Perp * OffsetDist);
        SpawnCylinderAt(MidPoint - Perp * OffsetDist);
    }
    else
    {
        // fallback: draw single cylinder for unknown order
        SpawnCylinderAt(MidPoint);
    }
}

// ----------------------------
// Element colors
// ----------------------------
FLinearColor APDBViewer::ElementColor(const FString& Element)
{
    if (Element.Equals("C",ESearchCase::IgnoreCase)) return FLinearColor(0.1f,0.1f,0.1f);
    if (Element.Equals("O",ESearchCase::IgnoreCase)) return FLinearColor::Red;
    if (Element.Equals("H",ESearchCase::IgnoreCase)||Element.Equals("D",ESearchCase::IgnoreCase)) return FLinearColor::White;
    if (Element.Equals("N",ESearchCase::IgnoreCase)) return FLinearColor::Blue;
    if (Element.Equals("S",ESearchCase::IgnoreCase)) return FLinearColor::Yellow;
    if (Element.Equals("CL",ESearchCase::IgnoreCase)) return FLinearColor(0.0f,1.0f,0.0f);
    if (Element.Equals("P",ESearchCase::IgnoreCase)) return FLinearColor(1.0f,0.5f,0.0f);
    if (Element.Equals("F",ESearchCase::IgnoreCase)) return FLinearColor(0.0f,1.0f,0.0f);
    if (Element.Equals("BR",ESearchCase::IgnoreCase)) return FLinearColor(0.6f,0.2f,0.2f);
    if (Element.Equals("I",ESearchCase::IgnoreCase)) return FLinearColor(0.4f,0.0f,0.8f);
    if (Element.Equals("FE",ESearchCase::IgnoreCase)) return FLinearColor(0.8f,0.4f,0.0f);
    if (Element.Equals("MG",ESearchCase::IgnoreCase)) return FLinearColor(0.0f,0.8f,0.0f);
    if (Element.Equals("ZN",ESearchCase::IgnoreCase)) return FLinearColor(0.5f,0.5f,0.5f);
    if (Element.Equals("CA",ESearchCase::IgnoreCase)) return FLinearColor(0.2f,0.6f,1.0f);
    if (Element.Equals("NA",ESearchCase::IgnoreCase)) return FLinearColor(0.0f,0.0f,1.0f);
    if (Element.Equals("K",ESearchCase::IgnoreCase)) return FLinearColor(0.5f,0.0f,1.0f);
    if (Element.Equals("CU",ESearchCase::IgnoreCase)) return FLinearColor(1.0f,0.5f,0.0f);
    if (Element.Equals("B",ESearchCase::IgnoreCase)) return FLinearColor(1.0f,0.7f,0.7f);
    return FLinearColor::Gray;
}



void APDBViewer::ClearCurrentStructure()
{
    // Destroy all atom meshes
    for (UStaticMeshComponent* Mesh : AllAtomMeshes)
    {
        if (Mesh && IsValid(Mesh))
        {
            Mesh->DestroyComponent();
        }
    }
    AllAtomMeshes.Empty();

    // Destroy all bond meshes
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

void APDBViewer::SaveStructureToFile(const FString& FilePath)
{
    if (CurrentPDBContent.IsEmpty())
    {
        UE_LOG(LogTemp, Warning, TEXT("No structure data to save"));
        return;
    }

    if (FFileHelper::SaveStringToFile(CurrentPDBContent, *FilePath))
    {
        UE_LOG(LogTemp, Log, TEXT("Structure saved to: %s"), *FilePath);
    }
    else
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to save structure to: %s"), *FilePath);
    }
}

void APDBViewer::LoadStructureFromFile(const FString& FilePath)
{
    FString FileContent;
    if (FFileHelper::LoadFileToString(FileContent, *FilePath))
    {
        UE_LOG(LogTemp, Log, TEXT("Loading structure from: %s"), *FilePath);
        
        // Clear existing structure
        ClearCurrentStructure();
        
        // Determine file type by extension
        FString Extension = FPaths::GetExtension(FilePath).ToLower();
        
        if (Extension == TEXT("pdb"))
        {
            CurrentPDBContent = FileContent;
            CurrentStructureID = FPaths::GetBaseFilename(FilePath);
            ParsePDB(FileContent);
        }
        else if (Extension == TEXT("cif"))
        {
            CurrentPDBContent = FileContent;
            CurrentStructureID = FPaths::GetBaseFilename(FilePath);
            ParseMMCIF(FileContent);
        }
        else
        {
            UE_LOG(LogTemp, Error, TEXT("Unsupported file format: %s"), *Extension);
        }
    }
    else
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to load file: %s"), *FilePath);
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
            OutFiles
        );
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
            OutFiles
        );
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


