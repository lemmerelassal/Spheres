#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Http.h"
#include "PDBViewer.generated.h"

// Structure to hold residue metadata during parsing
USTRUCT()
struct FResidueMetadata
{
    GENERATED_BODY()
    
    FString ResidueName;
    FString ResidueSeq;
    FString Chain;
    FString RecordType;
    
    FResidueMetadata()
        : ResidueName(TEXT(""))
        , ResidueSeq(TEXT(""))
        , Chain(TEXT(""))
        , RecordType(TEXT(""))
    {}
};

// Structure to hold residue information
USTRUCT()
struct FResidueInfo
{
    GENERATED_BODY()
    
    FString RecordType;   // e.g., "HETATM"
    FString Chain;        // e.g., "A"
    FString ResidueName;  // e.g., "SER"
    FString ResidueSeq;   // e.g., "1313"
    TArray<UStaticMeshComponent*> AtomMeshes;
    TArray<UStaticMeshComponent*> BondMeshes;
    bool bIsVisible;
    
    FResidueInfo()
        : RecordType(TEXT(""))
        , Chain(TEXT(""))
        , ResidueName(TEXT(""))
        , ResidueSeq(TEXT(""))
        , bIsVisible(true)
    {}
};

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
    
    // Batch residue creation
    void CreateResiduesFromAtomData(
        const TMap<FString, TMap<FString, FVector>>& ResidueAtoms,
        const TMap<FString, FResidueMetadata>& Metadata);
    
    // Bond parsing
    void FetchLigandBondsForResidue(
        const FString& ResidueKey,
        const FString& ResidueName,
        const TMap<FString, FVector>& AtomPositions);
    
    void ParseLigandCIFForResidue(
        const FString& FileContent,
        const TMap<FString, FVector>& AtomPositions,
        FResidueInfo* ResInfo);
    
    // Helper functions
    FString NormalizeAtomID(const FString& In) const;
    int32 ParseBondOrder(const FString& OrderStr) const;

    // Drawing
    void DrawSphere(
        float X, float Y, float Z,
        const FLinearColor& Color,
        USceneComponent* Parent,
        TArray<UStaticMeshComponent*>& OutArray);
    
    void DrawBond(
        const FVector& Start,
        const FVector& End,
        int32 Order,
        const FLinearColor& Color,
        USceneComponent* Parent,
        TArray<UStaticMeshComponent*>& OutArray);
    
    FLinearColor ElementColor(const FString& Element) const;

    // Storage for spawned components
    UPROPERTY()
    TArray<UStaticMeshComponent*> AllAtomMeshes;
    
    UPROPERTY()
    TArray<UStaticMeshComponent*> AllBondMeshes;

    // Current structure data
    FString CurrentPDBContent;
    FString CurrentStructureID;
    
    // Map from residue key (e.g., "SER_1313_A") to residue info
    TMap<FString, FResidueInfo*> ResidueMap;
    
    // Cleanup helper
    void ClearResidueMap();

public:
    // Event delegate
    DECLARE_DYNAMIC_MULTICAST_DELEGATE(FOnResiduesLoaded);

    UPROPERTY(BlueprintAssignable, Category = "PDB Viewer")
    FOnResiduesLoaded OnResiduesLoaded;

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
    
    // Residue visibility functions
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void ToggleResidueVisibility(const FString& ResidueKey);
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    TArray<FString> GetResidueList() const;
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    FString GetResidueDisplayName(const FString& ResidueKey) const;

    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    TArray<FString> GetLigandList() const;

private:
    // Helper for file dialogs
    bool ShowFileDialog(bool bSave, FString& OutFilePath);
};