#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "BlueSphere.generated.h"

USTRUCT()
struct FResidueData
{
    GENERATED_BODY()

    UPROPERTY()
    FString ResidueName;

    UPROPERTY()
    TArray<UStaticMeshComponent*> AtomSpheres;

    UPROPERTY()
    TArray<UStaticMeshComponent*> BondCylinders;
};

UCLASS()
class SPHERES_API ABlueSphere : public AActor
{
    GENERATED_BODY()

public:
    ABlueSphere();

protected:
    virtual void BeginPlay() override;

public:
    virtual void Tick(float DeltaTime) override;

    void LoadMoleculeFromJSON(const FString& FilePath);

    // Show/hide a specific residue
    void ToggleResidueVisibility(int32 ResidueIndex, bool bVisible);

    // Getter
    const TArray<FResidueData>& GetResidues() const { return Residues; }

private:
    void DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    void DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    FLinearColor ElementColor(const FString& Element);

private:
    UPROPERTY()
    UStaticMesh* SphereMeshAsset;

    UPROPERTY()
    UStaticMesh* CylinderMeshAsset;

    UPROPERTY()
    UMaterial* SphereMaterialAsset;

    UPROPERTY()
    USceneComponent* RootScene;

    // All residues
    UPROPERTY()
    TArray<FResidueData> Residues;
};
