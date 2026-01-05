#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Http.h"
#include "PDBViewer.generated.h"

USTRUCT()
struct FResidueMetadata
{
    GENERATED_BODY()
    FString ResidueName, ResidueSeq, Chain, RecordType;
};

USTRUCT()
struct FResidueInfo
{
    GENERATED_BODY()
    FString RecordType, Chain, ResidueName, ResidueSeq;
    TArray<UStaticMeshComponent*> AtomMeshes, BondMeshes;
    bool bIsVisible = true;
};

UCLASS()
class SPHERES_API APDBViewer : public AActor
{
    GENERATED_BODY()

public:
    APDBViewer();
    
    DECLARE_DYNAMIC_MULTICAST_DELEGATE(FOnResiduesLoaded);
    UPROPERTY(BlueprintAssignable, Category = "PDB Viewer") FOnResiduesLoaded OnResiduesLoaded;
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void SaveStructureToFile(const FString& FilePath);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void LoadStructureFromFile(const FString& FilePath);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void ClearCurrentStructure();
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void OpenSaveDialog();
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void OpenLoadDialog();
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void ToggleResidueVisibility(const FString& ResidueKey);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") TArray<FString> GetResidueList() const;
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") FString GetResidueDisplayName(const FString& ResidueKey) const;
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") TArray<FString> GetLigandList() const;

protected:
    virtual void BeginPlay() override;
    
    UPROPERTY() UStaticMesh* SphereMeshAsset;
    UPROPERTY() UStaticMesh* CylinderMeshAsset;
    UPROPERTY() UMaterial* SphereMaterialAsset;
    UPROPERTY() TArray<UStaticMeshComponent*> AllAtomMeshes;
    UPROPERTY() TArray<UStaticMeshComponent*> AllBondMeshes;
    
    FString CurrentPDBContent, CurrentStructureID;
    TMap<FString, FResidueInfo*> ResidueMap;
    
    void FetchAndDisplayStructure(const FString& PDB_ID);
    void FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> Callback);
    void ParsePDB(const FString& FileContent);
    void ParseMMCIF(const FString& FileContent);
    void CreateResiduesFromAtomData(const TMap<FString, TMap<FString, FVector>>& ResidueAtoms, const TMap<FString, FResidueMetadata>& Metadata);
    void FetchLigandBondsForResidue(const FString& ResidueKey, const FString& ResidueName, const TMap<FString, FVector>& AtomPositions);
    void ParseLigandCIFForResidue(const FString& FileContent, const TMap<FString, FVector>& AtomPositions, FResidueInfo* ResInfo);
    FString NormalizeAtomID(const FString& In) const;
    int32 ParseBondOrder(const FString& OrderStr) const;
    void DrawSphere(float X, float Y, float Z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    void DrawBond(const FVector& Start, const FVector& End, int32 Order, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    FLinearColor GetElementColor(const FString& Element) const;
    void ClearResidueMap();
    bool ShowFileDialog(bool bSave, FString& OutFilePath);
};