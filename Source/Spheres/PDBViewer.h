#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Http.h"
#include "PDBViewer.generated.h"

UCLASS()
class SPHERES_API APDBViewer : public AActor
{
    GENERATED_BODY()

public:
    APDBViewer();

protected:
    virtual void BeginPlay() override;

    // Assets
    UPROPERTY()
    UStaticMesh* SphereMeshAsset;
    
    UPROPERTY()
    UStaticMesh* CylinderMeshAsset;
    
    UPROPERTY()
    UMaterial* SphereMaterialAsset;

    // Parsing & Fetching
    void FetchAndDisplayStructure(const FString& PDB_ID);
    void FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> Callback);
    void ParsePDB(const FString& FileContent);
    void ParseMMCIF(const FString& FileContent);
    void FetchLigandCIF(const FString& ResidueName, const TMap<FString, FVector>& AtomPositions);
    void ParseLigandCIF(const FString& FileContent, const TMap<FString, FVector>& AtomPositions);

    // Drawing
    void DrawSphere(float x, float y, float z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    void DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    FLinearColor ElementColor(const FString& Element);

    // Storage for spawned components
    UPROPERTY()
    TArray<UStaticMeshComponent*> AllAtomMeshes;
    
    UPROPERTY()
    TArray<UStaticMeshComponent*> AllBondMeshes;

    // Current structure data
    FString CurrentPDBContent;
    FString CurrentStructureID;

public:
    // Save/Load functions
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void SaveStructureToFile(const FString& FilePath);
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void LoadStructureFromFile(const FString& FilePath);
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void ClearCurrentStructure();
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void OpenSaveDialog();
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void OpenLoadDialog();

private:
    // Helper for file dialogs
    bool ShowFileDialog(bool bSave, FString& OutFilePath);
};