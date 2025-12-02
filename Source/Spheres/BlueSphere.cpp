
// BlueSphere.cpp

#include "DrawDebugHelpers.h"  // Needed for DrawDebugLine

#include "BlueSphere.h"
#include "Components/StaticMeshComponent.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "UObject/ConstructorHelpers.h"
#include "Misc/FileHelper.h"
#include "Misc/Paths.h"
#include "Serialization/JsonReader.h"
#include "Serialization/JsonSerializer.h"

ABlueSphere::ABlueSphere()
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

void ABlueSphere::BeginPlay()
{
    Super::BeginPlay();

    SetActorLocation(FVector::ZeroVector);

    FString JSONPath = TEXT("D:/Golang/pdbParserFinal/5ENB.json");
    LoadMoleculeFromJSON(JSONPath);
}

void ABlueSphere::LoadMoleculeFromJSON(const FString& FilePath)
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

            // Residue level
            for (auto& ResPair : ChainObj->Values)
            {
                TSharedPtr<FJsonObject> ResObj = ResPair.Value->AsObject();
                if (!ResObj.IsValid()) continue;

                const TArray<TSharedPtr<FJsonValue>>* AtomsArray;
                if (!ResObj->TryGetArrayField("atoms", AtomsArray)) continue;

                // Store atom positions for bonds
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

                    DrawSphere(Pos.X, Pos.Y, Pos.Z, ElementColor(Element), RootComponent);
                }

                // Bonds
                const TArray<TSharedPtr<FJsonValue>>* BondsArray;
                if (ResObj->TryGetArrayField("bonds", BondsArray))
                {
                    for (const TSharedPtr<FJsonValue>& BondValue : *BondsArray)
                    {
                        TSharedPtr<FJsonObject> BondObj = BondValue->AsObject();
                        if (!BondObj.IsValid()) continue;

                        int32 From = BondObj->GetIntegerField("from");
                        int32 To = BondObj->GetIntegerField("to");
                        int32 Order = BondObj->GetIntegerField("order");

                        if(Order > 1) {
                            UE_LOG(LogTemp, Warning, TEXT("Multiple bond orders not supported yet."));
                        }

                        if (AtomPositions.IsValidIndex(From) && AtomPositions.IsValidIndex(To))
                        {
                            DrawBond(AtomPositions[From], AtomPositions[To], Order, FLinearColor::Gray, RootComponent);
                        }
                    }
                }
            }
        }
    }
}

FLinearColor ABlueSphere::ElementColor(const FString& Element)
{
    if (Element == "C") return FLinearColor(0.1f, 0.1f, 0.1f);
    if (Element == "O") return FLinearColor::Red;
    if (Element == "H") return FLinearColor::White;
    if (Element == "N") return FLinearColor::Blue;
    return FLinearColor::Green;
}

void ABlueSphere::DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent)
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
}



void ABlueSphere::DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent)
{
    FVector BondVector = End - Start;
    float Length = BondVector.Size();
    FVector MidPoint = Start + 0.5f * BondVector;
    FRotator Rotation = FRotationMatrix::MakeFromZ(BondVector).Rotator();

    float DefaultHalfHeight = 50.0f;
    float ZScale = Length / (2.0f * DefaultHalfHeight);
    float OffsetAmount = 8.0f;
    FVector UpVector = FVector::UpVector;
    if (FMath::Abs(FVector::DotProduct(BondVector.GetSafeNormal(), UpVector)) > 0.99f)
        UpVector = FVector::RightVector;
    FVector OffsetDir = FVector::CrossProduct(BondVector, UpVector).GetSafeNormal();

    for (int32 i = 0; i < Order; i++)
    {
        FVector Offset = FVector::ZeroVector;
        if (Order == 2)
            Offset = OffsetDir * ((i == 0) ? OffsetAmount : -OffsetAmount);
        else if (Order == 3)
            Offset = OffsetDir * ((i - 1) * OffsetAmount);

        UStaticMeshComponent* Cylinder = NewObject<UStaticMeshComponent>(this);
        Cylinder->SetStaticMesh(CylinderMeshAsset);
        Cylinder->SetWorldLocation(MidPoint + Offset);
        Cylinder->SetWorldRotation(Rotation);
        Cylinder->SetWorldScale3D(FVector(0.1f, 0.1f, ZScale));

        Cylinder->RegisterComponent();
        AddInstanceComponent(Cylinder);
        Cylinder->AttachToComponent(Parent, FAttachmentTransformRules::KeepWorldTransform);

        UMaterialInstanceDynamic* Mat = UMaterialInstanceDynamic::Create(SphereMaterialAsset, this);
        Mat->SetVectorParameterValue("Color", Color);
        Cylinder->SetMaterial(0, Mat);
    }
}



void ABlueSphere::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);
}
