#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Ligands.generated.h"

USTRUCT()
struct FLigandData
{
    GENERATED_BODY()

    UPROPERTY()
    FString LigandName;

    UPROPERTY()
    TArray<UStaticMeshComponent*> AtomSpheres;

    UPROPERTY()
    TArray<UStaticMeshComponent*> BondCylinders;
};

UCLASS()
class SPHERES_API ALigands : public AActor
{
    GENERATED_BODY()

public:
    ALigands();

protected:
    virtual void BeginPlay() override;

public:
    virtual void Tick(float DeltaTime) override;

    // Load ligand data from a JSON file
    void LoadMoleculeFromJSON(const FString& FilePath);

    // Show or hide a specific ligand
    void ToggleLigandVisibility(int32 LigandIndex, bool bVisible);

    // Get reference to all ligands
    const TArray<FLigandData>& GetLigands() const;   // <-- declaration only (no body here)

private:
    // Drawing helpers
    void DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    void DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    FLinearColor ElementColor(const FString& Element);

private:
    UPROPERTY()
    USceneComponent* RootScene;

    UPROPERTY()
    UStaticMesh* SphereMeshAsset;

    UPROPERTY()
    UStaticMesh* CylinderMeshAsset;

    UPROPERTY()
    UMaterial* SphereMaterialAsset;

    // All ligands
    UPROPERTY()
    TArray<FLigandData> LigandsArray;
};
