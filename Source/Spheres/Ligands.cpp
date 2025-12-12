#include "Ligands.h"
#include "LigandVisibilityWidget.h"
#include "Blueprint/WidgetBlueprintLibrary.h"
#include "DrawDebugHelpers.h"
#include "Components/StaticMeshComponent.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "UObject/ConstructorHelpers.h"
#include "Misc/FileHelper.h"
#include "Misc/Paths.h"
#include "Serialization/JsonReader.h"
#include "Serialization/JsonSerializer.h"

ALigands::ALigands()
{
    PrimaryActorTick.bCanEverTick = true;

    // Sphere mesh
    static ConstructorHelpers::FObjectFinder<UStaticMesh> SphereRef(TEXT("/Engine/BasicShapes/Sphere.Sphere"));
    if (SphereRef.Succeeded()) SphereMeshAsset = SphereRef.Object;

    // Cylinder mesh (for bonds)
    static ConstructorHelpers::FObjectFinder<UStaticMesh> CylinderRef(TEXT("/Engine/BasicShapes/Cylinder.Cylinder"));
    if (CylinderRef.Succeeded()) CylinderMeshAsset = CylinderRef.Object;

    // Material
    static ConstructorHelpers::FObjectFinder<UMaterial> MatRef(TEXT("/Engine/BasicShapes/BasicShapeMaterial.BasicShapeMaterial"));
    if (MatRef.Succeeded()) SphereMaterialAsset = MatRef.Object;

    // Root
    RootScene = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
    RootComponent = RootScene;
}

void ALigands::BeginPlay()
{
    Super::BeginPlay();

    SetActorLocation(FVector::ZeroVector);

    FString JSONPath = TEXT("D:/Golang/attach_bonds/output_with_bonds.json");
    LoadMoleculeFromJSON(JSONPath);

    // Create and display the widget
    TSharedPtr<FLigandVisibilityWidget> LigandWidget = SNew(FLigandVisibilityWidget)
        .Ligands(this);

    GEngine->GameViewport->AddViewportWidgetContent(LigandWidget.ToSharedRef());
}

void ALigands::LoadMoleculeFromJSON(const FString& FilePath)
{
    FString JsonText;
    if (!FFileHelper::LoadFileToString(JsonText, *FilePath))
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to load JSON: %s"), *FilePath);
        return;
    }

    TSharedPtr<FJsonValue> RootValue;
    TSharedRef<TJsonReader<>> Reader = TJsonReaderFactory<>::Create(JsonText);

    if (!FJsonSerializer::Deserialize(Reader, RootValue) || !RootValue.IsValid())
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to parse JSON"));
        return;
    }

    const float Scale = 50.f;

    // Helper to parse a single ligand JSON object
    auto ParseSingleLigand = [&](TSharedPtr<FJsonObject> LigObj)
    {
        if (!LigObj.IsValid()) return; // nothing to do

        FLigandData NewLigand;

        // Build a friendly ligand name if fields exist
        FString PDB = LigObj->HasField(TEXT("pdb")) ? LigObj->GetStringField(TEXT("pdb")) : FString();
        FString Chain = LigObj->HasField(TEXT("chain")) ? LigObj->GetStringField(TEXT("chain")) : FString();
        FString LigandID = LigObj->HasField(TEXT("ligand")) ? LigObj->GetStringField(TEXT("ligand")) : FString(TEXT("Ligand"));

        if (!PDB.IsEmpty() && !Chain.IsEmpty())
        {
            NewLigand.LigandName = FString::Printf(TEXT("%s (%s Chain %s)"), *LigandID, *PDB, *Chain);
        }
        else if (!PDB.IsEmpty())
        {
            NewLigand.LigandName = FString::Printf(TEXT("%s (%s)"), *LigandID, *PDB);
        }
        else
        {
            NewLigand.LigandName = LigandID;
        }

        // Get atoms object
        if (!LigObj->HasField(TEXT("atoms")))
        {
            UE_LOG(LogTemp, Warning, TEXT("Ligand %s missing \"atoms\" field, skipping"), *NewLigand.LigandName);
            return;
        }

        TSharedPtr<FJsonObject> AtomsObject = LigObj->GetObjectField(TEXT("atoms"));
        if (!AtomsObject.IsValid())
        {
            UE_LOG(LogTemp, Warning, TEXT("Ligand %s has invalid \"atoms\" object, skipping"), *NewLigand.LigandName);
            return;
        }

        TArray<FString> AtomKeys;
        AtomsObject->Values.GetKeys(AtomKeys);

        // Skip ligands with fewer than 6 atoms
        if (AtomKeys.Num() < 6)
        {
            UE_LOG(LogTemp, Warning, TEXT("Skipping ligand %s (only %d atoms)"), *NewLigand.LigandName, AtomKeys.Num());
            return;
        }

        // Map atom name -> index (populate after validating count)
        TMap<FString, int32> NameToIndex;
        NameToIndex.Reserve(AtomKeys.Num());

        // Now create spheres and populate atom arrays
        for (const FString& AtomName : AtomKeys)
        {
            TSharedPtr<FJsonObject> AtomObj = AtomsObject->GetObjectField(AtomName);
            if (!AtomObj.IsValid()) continue;

            double X = 0.0, Y = 0.0, Z = 0.0;
            if (AtomObj->HasField(TEXT("x"))) X = AtomObj->GetNumberField(TEXT("x"));
            if (AtomObj->HasField(TEXT("y"))) Y = AtomObj->GetNumberField(TEXT("y"));
            if (AtomObj->HasField(TEXT("z"))) Z = AtomObj->GetNumberField(TEXT("z"));

            FString Element = AtomObj->HasField(TEXT("element")) ? AtomObj->GetStringField(TEXT("element")) : FString(TEXT("C"));

            FVector Pos = FVector(X * Scale, Y * Scale, Z * Scale);

            DrawSphere(Pos.X, Pos.Y, Pos.Z, ElementColor(Element), RootComponent, NewLigand.AtomSpheres);

            NewLigand.AtomPositions.Add(Pos);
            NewLigand.AtomElements.Add(Element);
            NewLigand.AtomNames.Add(AtomName);

            NameToIndex.Add(AtomName, NewLigand.AtomPositions.Num() - 1);
        }

        // Parse bonds array: each bond is expected to be an object with "atom1","atom2", and optional "order"
        const TArray<TSharedPtr<FJsonValue>>* BondsArray = nullptr;
        if (LigObj->TryGetArrayField(TEXT("bonds"), BondsArray) && BondsArray)
        {
            for (const TSharedPtr<FJsonValue>& BondVal : *BondsArray)
            {
                if (!BondVal.IsValid()) continue;

                // allow object form { "atom1":"A", "atom2":"B", "order":1 } OR array form ["A","B",1]
                if (BondVal->Type == EJson::Object)
                {
                    TSharedPtr<FJsonObject> BondObj = BondVal->AsObject();
                    if (!BondObj.IsValid()) continue;

                    FString AName = BondObj->HasField(TEXT("atom1")) ? BondObj->GetStringField(TEXT("atom1")) : FString();
                    FString BName = BondObj->HasField(TEXT("atom2")) ? BondObj->GetStringField(TEXT("atom2")) : FString();
                    int32 Order = BondObj->HasField(TEXT("order")) ? BondObj->GetIntegerField(TEXT("order")) : 1;

                    if (AName.IsEmpty() || BName.IsEmpty()) continue;

                    int32 AIndex = NameToIndex.Contains(AName) ? NameToIndex[AName] : INDEX_NONE;
                    int32 BIndex = NameToIndex.Contains(BName) ? NameToIndex[BName] : INDEX_NONE;

                    if (AIndex != INDEX_NONE && BIndex != INDEX_NONE)
                    {
                        DrawBond(NewLigand.AtomPositions[AIndex], NewLigand.AtomPositions[BIndex], Order, FLinearColor::Gray, RootComponent, NewLigand.BondCylinders);
                    }
                    else
                    {
                        UE_LOG(LogTemp, Warning, TEXT("Bond references missing atom(s): %s - %s in ligand %s"), *AName, *BName, *NewLigand.LigandName);
                    }
                }
                else if (BondVal->Type == EJson::Array)
                {
                    const TArray<TSharedPtr<FJsonValue>>& Pair = BondVal->AsArray();
                    if (Pair.Num() < 2) continue;

                    FString AName = Pair[0]->AsString();
                    FString BName = Pair[1]->AsString();
                    int32 Order = (Pair.Num() >= 3) ? (int32)Pair[2]->AsNumber() : 1;

                    int32 AIndex = NameToIndex.Contains(AName) ? NameToIndex[AName] : INDEX_NONE;
                    int32 BIndex = NameToIndex.Contains(BName) ? NameToIndex[BName] : INDEX_NONE;

                    if (AIndex != INDEX_NONE && BIndex != INDEX_NONE)
                    {
                        DrawBond(NewLigand.AtomPositions[AIndex], NewLigand.AtomPositions[BIndex], Order, FLinearColor::Gray, RootComponent, NewLigand.BondCylinders);
                    }
                    else
                    {
                        UE_LOG(LogTemp, Warning, TEXT("Bond references missing atom(s): %s - %s in ligand %s"), *AName, *BName, *NewLigand.LigandName);
                    }
                }
            }
        }

        // Store the ligand
        LigandsArray.Add(NewLigand);
    };

    // If root is an array of ligands, parse each entry
    if (RootValue->Type == EJson::Array)
    {
        const TArray<TSharedPtr<FJsonValue>>& RootArray = RootValue->AsArray();
        for (const TSharedPtr<FJsonValue>& Element : RootArray)
        {
            if (!Element.IsValid()) continue;
            if (Element->Type == EJson::Object)
            {
                ParseSingleLigand(Element->AsObject());
            }
        }
        return;
    }

    // If root is an object, treat it as one ligand
    if (RootValue->Type == EJson::Object)
    {
        ParseSingleLigand(RootValue->AsObject());
        return;
    }

    UE_LOG(LogTemp, Error, TEXT("Unsupported JSON root type for molecule file"));
}

FLinearColor ALigands::ElementColor(const FString& Element)
{
    if (Element.IsEmpty()) return FLinearColor::Black;


    if(Element == "C") return FLinearColor(0.1f, 0.1f, 0.1f);
    if(Element == "O") return FLinearColor::Red;
    if(Element == "H") return FLinearColor::White;
    if(Element == "D") return FLinearColor::White;
    if(Element == "N") return FLinearColor::Blue;
    if(Element == "S") return FLinearColor::Yellow;
    if(Element == "CL") return FLinearColor(0.0f, 1.0f, 0.0f);
    if(Element == "P") return FLinearColor(1.0f, 0.5f, 0.0f);
    if(Element == "F") return FLinearColor(0.0f, 1.0f, 0.0f);
    if(Element == "BR") return FLinearColor(0.6f, 0.2f, 0.2f);
    if(Element == "I") return FLinearColor(0.4f, 0.0f, 0.8f);
    if(Element == "FE") return FLinearColor(0.8f, 0.4f, 0.0f);
    if(Element == "MG") return FLinearColor(0.0f, 0.8f, 0.0f);
    if(Element == "ZN") return FLinearColor(0.5f, 0.5f, 0.5f);
    if(Element == "CA") return FLinearColor(0.2f, 0.6f, 1.0f);
    if(Element == "NA") return FLinearColor(0.0f, 0.0f, 1.0f);
    if(Element == "K") return FLinearColor(0.5f, 0.0f, 1.0f);
    if(Element == "CU") return FLinearColor(1.0f, 0.5f, 0.0f);
    if(Element == "B") return FLinearColor(1.0f, 0.7f, 0.7f);
    return FLinearColor::Gray;
    

}

void ALigands::DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray)
{
    // Main sphere (colored, slightly flat/emissive to approximate cel shading)
    UStaticMeshComponent* Sphere = NewObject<UStaticMeshComponent>(this);
    Sphere->RegisterComponent();
    Sphere->AttachToComponent(Parent, FAttachmentTransformRules::KeepRelativeTransform);
    Sphere->SetStaticMesh(SphereMeshAsset);
    Sphere->SetRelativeLocation(FVector(x, y, z));
    Sphere->SetWorldScale3D(FVector(0.5f));
    Sphere->SetCollisionEnabled(ECollisionEnabled::NoCollision);
    Sphere->SetCastShadow(true);

    UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
    Mat->SetVectorParameterValue(FName("Color"), Color);

    // Boost emissive a bit so shading looks flatter (depends on the base material)
    Mat->SetVectorParameterValue(FName("EmissiveColor"), Color * 0.35f);
    Mat->SetScalarParameterValue(FName("EmissiveIntensity"), 5.0f);

    // If your base material supports parameters like roughness/specular, you can
    // set them here to reduce specular highlights and favor flat diffuse.
    // Mat->SetScalarParameterValue(FName("Roughness"), 1.0f);
    // Mat->SetScalarParameterValue(FName("Specular"), 0.0f);

    Sphere->SetMaterial(0, Mat);

    // Optional: enable custom depth if you later want a post-process outline
    Sphere->SetRenderCustomDepth(false);

    OutArray.Add(Sphere);

    // Outline/backface sphere (slightly larger, reverse-cull so backfaces render)
    UStaticMeshComponent* Outline = NewObject<UStaticMeshComponent>(this);
    Outline->RegisterComponent();
    Outline->AttachToComponent(Parent, FAttachmentTransformRules::KeepRelativeTransform);
    Outline->SetStaticMesh(SphereMeshAsset);
    Outline->SetRelativeLocation(FVector(x, y, z));

    // Slightly larger to produce a silhouette rim
    const float OutlineScale = 0.54f;
    Outline->SetWorldScale3D(FVector(OutlineScale));
    Outline->SetCollisionEnabled(ECollisionEnabled::NoCollision);
    Outline->SetCastShadow(false);

    // Invert culling so the backfaces of the scaled-up sphere show as a rim
    Outline->bReverseCulling = true;

    UMaterialInstanceDynamic* OutlineMat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
    OutlineMat->SetVectorParameterValue(FName("Color"), FLinearColor::Black);
    // Make it emissive so it's always solid black regardless of lighting
    OutlineMat->SetVectorParameterValue(FName("EmissiveColor"), FLinearColor::Black);
    OutlineMat->SetScalarParameterValue(FName("EmissiveIntensity"), 50.0f);

    Outline->SetMaterial(0, OutlineMat);

    OutArray.Add(Outline);
}

void ALigands::DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray)
{
    FVector BondVector = End - Start;
    float Length = BondVector.Size();
    FVector MidPoint = Start + 0.5f * BondVector;
    FRotator Rotation = FRotationMatrix::MakeFromZ(BondVector).Rotator();

    const float DefaultHalfHeight = 50.0f;
    const float ZScale = Length / (2.0f * DefaultHalfHeight);

    UStaticMeshComponent* Cylinder = NewObject<UStaticMeshComponent>(this);
    Cylinder->SetStaticMesh(CylinderMeshAsset);
    Cylinder->SetWorldLocation(MidPoint);
    Cylinder->SetWorldRotation(Rotation);
    Cylinder->SetWorldScale3D(FVector(0.1f, 0.1f, ZScale));
    Cylinder->RegisterComponent();
    Cylinder->AttachToComponent(Parent, FAttachmentTransformRules::KeepWorldTransform);

    UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
    Mat->SetVectorParameterValue("Color", Color);
    Cylinder->SetMaterial(0, Mat);

    OutArray.Add(Cylinder);
}

void ALigands::ToggleLigandVisibility(int32 LigandIndex, bool bVisible)
{
    if (!LigandsArray.IsValidIndex(LigandIndex)) return;

    FLigandData& Ligand = LigandsArray[LigandIndex];

    for (UStaticMeshComponent* Atom : Ligand.AtomSpheres)
    {
        if (Atom) Atom->SetVisibility(bVisible);
    }

    for (UStaticMeshComponent* Bond : Ligand.BondCylinders)
    {
        if (Bond) Bond->SetVisibility(bVisible);
    }
}

void ALigands::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);
}

const TArray<FLigandData>& ALigands::GetLigands() const
{
    return LigandsArray;
}
