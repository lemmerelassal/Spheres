#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "BlueSphere.generated.h"

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

    
    // Methods to show/hide residues
    void ToggleResidueVisibility(bool bVisible);

    
    // Getter for AtomSpheres
    TArray<UStaticMeshComponent*>& GetAtomSpheres() { return AtomSpheres; }



private:
    // Sphere
    void DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent);
    void DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent);

    FLinearColor ElementColor(const FString& Element);

    UPROPERTY()
    UStaticMesh* SphereMeshAsset;

    UPROPERTY()
    UStaticMesh* CylinderMeshAsset;

    UPROPERTY()
    UMaterial* SphereMaterialAsset;

    UPROPERTY()
    USceneComponent* RootScene;

    
    TArray<UStaticMeshComponent*> AtomSpheres;  // To hold the residue sphere components
    TArray<UStaticMeshComponent*> BondCylinders; // To hold the bond cylinder components
};
