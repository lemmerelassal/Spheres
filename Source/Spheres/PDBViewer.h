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

    UFUNCTION(BlueprintCallable, Category="Molecule")
    void FetchAndDisplayStructure(const FString& PDB_ID);

protected:
    virtual void BeginPlay() override;

private:
    /* ---------- NETWORKING ---------- */
    void FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> Callback);

    /* ---------- MAIN STRUCTURE PARSING ---------- */
    void ParsePDB(const FString& FileContent);
    void ParseMMCIF(const FString& FileContent);

    /* ---------- LIGAND TOPOLOGY (CIF) ---------- */
    void FetchLigandCIF(const FString& ResidueName,
                        const TMap<FString, FVector>& AtomPositions,
                        const FVector& Offset);

    void ParseLigandCIF(const FString& FileContent,
                        const TMap<FString, FVector>& AtomPositions,
                        const FVector& Offset);

    /* ---------- DRAWING ---------- */
    FLinearColor ElementColor(const FString& Element);

    void DrawSphere(float x, float y, float z,
                    const FLinearColor& Color,
                    USceneComponent* Parent,
                    TArray<UStaticMeshComponent*>& OutArray);

    void DrawBond(const FVector& Start,
                  const FVector& End,
                  int32 Order,
                  const FLinearColor& Color,
                  USceneComponent* Parent,
                  TArray<UStaticMeshComponent*>& OutArray);

    /* ---------- ASSETS ---------- */
    // Assign these in the Unreal Editor or constructor.
    UPROPERTY(EditAnywhere, Category="Molecule|Meshes")
    UStaticMesh* SphereMeshAsset;

    UPROPERTY(EditAnywhere, Category="Molecule|Meshes")
    UStaticMesh* CylinderMeshAsset;

    UPROPERTY(EditAnywhere, Category="Molecule|Materials")
    UMaterialInterface* SphereMaterialAsset;
};
