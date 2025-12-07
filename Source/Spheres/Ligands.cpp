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

    FString JSONPath = TEXT("D:/Golang/pdb-group-parser/ligands.json");
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

        // Some ID level (like "5S8M")
        for (auto& IDPair : PDBObj->Values)
        {
            TSharedPtr<FJsonObject> IDObj = IDPair.Value->AsObject();
            if (!IDObj.IsValid()) continue;

            // Chain level (like "A")
            for (auto& ChainPair : IDObj->Values)
            {
                const TArray<TSharedPtr<FJsonValue>>* AtomArray;
                if (!ChainPair.Value->TryGetArray(AtomArray)) continue;

                FLigandData NewLigand;
                NewLigand.LigandName = FString::Printf(TEXT("%s %s (%s)"), *IDPair.Key, *ChainPair.Key, *PDBPair.Key);

                TArray<FVector> AtomPositions;
                TArray<FString> AtomElements;

                for (const TSharedPtr<FJsonValue>& AtomValue : *AtomArray)
                {
                    TSharedPtr<FJsonObject> AtomObj = AtomValue->AsObject();
                    if (!AtomObj.IsValid()) continue;

                    FString Element = AtomObj->GetStringField("hetatm_name");

                    TSharedPtr<FJsonObject> PosObj = AtomObj->GetObjectField("position");
                    double X = PosObj->GetNumberField("X");
                    double Y = PosObj->GetNumberField("Y");
                    double Z = PosObj->GetNumberField("Z");

                    FVector Pos = FVector(X * Scale, Y * Scale, Z * Scale);
                    AtomPositions.Add(Pos);
                    AtomElements.Add(Element);

                    DrawSphere(Pos.X, Pos.Y, Pos.Z, ElementColor(Element), RootComponent, NewLigand.AtomSpheres);
                }

                // Bonds handling could be added here if available

                LigandsArray.Add(NewLigand);
            }
        }
    }
}


FLinearColor ALigands::ElementColor(const FString& Element)
{
    if (Element.IsEmpty()) return FLinearColor::Black;

TCHAR FirstChar = Element[0];



    if(Element.StartsWith("CL")) return FLinearColor(0.0f, 1.0f, 0.0f);
    if(Element.StartsWith("NA")) return FLinearColor(0.0f, 0.0f, 1.0f);
    if(Element.StartsWith("MG")) return FLinearColor(0.0f, 0.8f, 0.0f);
    if(Element.StartsWith("CA")) return FLinearColor(0.5f, 0.5f, 0.5f);
    if(Element.StartsWith("FE")) return FLinearColor(0.8f, 0.4f, 0.0f);
    if(Element.StartsWith("BR")) return FLinearColor(0.6f, 0.2f, 0.2f);
    if(Element.StartsWith("F")) return FLinearColor(0.0f, 1.0f, 0.0f);
    if(Element.StartsWith("I")) return FLinearColor(0.4f, 0.0f, 0.8f);
    if(Element.StartsWith("K")) return FLinearColor(0.5f, 0.0f, 1.0f);
    if(Element.StartsWith("ZN")) return FLinearColor(0.5f, 0.5f, 0.0f);
    if(Element.StartsWith("CU")) return FLinearColor(0.8f, 0.5f, 0.2f);
    if(Element.StartsWith("SE")) return FLinearColor(1.0f, 0.5f, 0.0f);

    switch (FirstChar)
    {
    case 'C': return FLinearColor(0.1f, 0.1f, 0.1f);
    case 'O': return FLinearColor::Red;
    case 'H': return FLinearColor::White;
    case 'N': return FLinearColor::Blue;
    case 'S': return FLinearColor::Yellow;
    case 'P': return FLinearColor(1.f, 0.5f, 0.f);
    default:  return FLinearColor::Gray;
    }

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
