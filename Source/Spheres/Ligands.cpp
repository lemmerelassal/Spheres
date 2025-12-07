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

    FString JSONPath = TEXT("D:/Golang/pdbParserFinal/5ENB.json");
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

    TSharedPtr<FJsonObject> JsonObject;
    TSharedRef<TJsonReader<>> Reader = TJsonReaderFactory<>::Create(JsonText);

    if (!FJsonSerializer::Deserialize(Reader, JsonObject) || !JsonObject.IsValid())
    {
        UE_LOG(LogTemp, Error, TEXT("Failed to parse JSON"));
        return;
    }

    const float Scale = 50.f;

    // PDB ID level
    for (auto& PDBPair : JsonObject->Values)
    {
        TSharedPtr<FJsonObject> PDBObj = PDBPair.Value->AsObject();
        if (!PDBObj.IsValid()) continue;

        // Chain level
        for (auto& ChainPair : PDBObj->Values)
        {
            TSharedPtr<FJsonObject> ChainObj = ChainPair.Value->AsObject();
            if (!ChainObj.IsValid()) continue;

            // Ligand level
            for (auto& LigandPair : ChainObj->Values)
            {
                TSharedPtr<FJsonObject> LigandObj = LigandPair.Value->AsObject();
                if (!LigandObj.IsValid()) continue;

                FString LigandID = LigandPair.Key;
                FString LigandName = TEXT("Unknown");
                if (LigandObj->HasField("resName"))
                {
                    LigandName = LigandObj->GetStringField("resName");
                }

                FString DisplayName = FString::Printf(TEXT("%s %s (%s)"), *ChainPair.Key, *LigandID, *LigandName);

                FLigandData NewLigand;
                NewLigand.LigandName = DisplayName;

                const TArray<TSharedPtr<FJsonValue>>* AtomsArray;
                if (!LigandObj->TryGetArrayField("atoms", AtomsArray)) continue;

                TArray<FVector> AtomPositions;
                TArray<FString> AtomElements;

                for (const TSharedPtr<FJsonValue>& AtomValue : *AtomsArray)
                {
                    TSharedPtr<FJsonObject> AtomObj = AtomValue->AsObject();
                    if (!AtomObj.IsValid()) continue;

                    FString Element = AtomObj->GetStringField("element");
                    double X = AtomObj->GetNumberField("x");
                    double Y = AtomObj->GetNumberField("y");
                    double Z = AtomObj->GetNumberField("z");

                    FVector Pos = FVector(X * Scale, Y * Scale, Z * Scale);
                    AtomPositions.Add(Pos);
                    AtomElements.Add(Element);

                    DrawSphere(Pos.X, Pos.Y, Pos.Z, ElementColor(Element), RootComponent, NewLigand.AtomSpheres);
                }

                // Bonds
                const TArray<TSharedPtr<FJsonValue>>* BondsArray;
                if (LigandObj->TryGetArrayField("bonds", BondsArray))
                {
                    for (const TSharedPtr<FJsonValue>& BondValue : *BondsArray)
                    {
                        TSharedPtr<FJsonObject> BondObj = BondValue->AsObject();
                        if (!BondObj.IsValid()) continue;

                        int32 From = BondObj->GetIntegerField("from");
                        int32 To = BondObj->GetIntegerField("to");
                        int32 Order = BondObj->GetIntegerField("order");

                        if (AtomPositions.IsValidIndex(From) && AtomPositions.IsValidIndex(To))
                        {
                            DrawBond(AtomPositions[From], AtomPositions[To], Order, FLinearColor::Gray, RootComponent, NewLigand.BondCylinders);
                        }
                    }
                }

                LigandsArray.Add(NewLigand);
            }
        }
    }
}

FLinearColor ALigands::ElementColor(const FString& Element)
{
    if (Element == "C") return FLinearColor(0.1f, 0.1f, 0.1f);
    if (Element == "O") return FLinearColor::Red;
    if (Element == "H") return FLinearColor::White;
    if (Element == "N") return FLinearColor::Blue;
    return FLinearColor::Green;
}

void ALigands::DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray)
{
    UStaticMeshComponent* Sphere = NewObject<UStaticMeshComponent>(this);
    Sphere->RegisterComponent();
    Sphere->AttachToComponent(Parent, FAttachmentTransformRules::KeepRelativeTransform);
    Sphere->SetStaticMesh(SphereMeshAsset);
    Sphere->SetRelativeLocation(FVector(x, y, z));
    Sphere->SetWorldScale3D(FVector(0.5f));

    UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
    Mat->SetVectorParameterValue(FName("Color"), Color);
    Mat->SetScalarParameterValue(FName("EmissiveIntensity"), 5.f);
    Sphere->SetMaterial(0, Mat);

    OutArray.Add(Sphere);
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
