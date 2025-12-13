// PDBViewer.h

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "PDBViewer.generated.h"

UCLASS()
class SPHERES_API APDBViewer : public AActor
{
    GENERATED_BODY()

public:
    APDBViewer();

protected:
    virtual void BeginPlay() override;

private:
    // Meshes and materials
    UStaticMesh* SphereMeshAsset;
    UStaticMesh* CylinderMeshAsset;
    UMaterial* SphereMaterialAsset;

    // Fetch & parse
    void FetchAndDisplayStructure(const FString& PDB_ID);
    void FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> Callback);

    void ParsePDB(const FString& FileContent);
    void ParseMMCIF(const FString& FileContent);

    void FetchLigandCIF(const FString& ResidueName, const TMap<FString, FVector>& AtomPositions);
    void ParseLigandCIF(const FString& FileContent, const TMap<FString, FVector>& AtomPositions);

    // Draw functions
    void DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    void DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);

    // Helpers
    FLinearColor ElementColor(const FString& Element);
};
