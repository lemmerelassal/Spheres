#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Http.h"
#include "MMGBSA.h"
#include "PDBViewer.generated.h"

USTRUCT()
struct FResidueMetadata
{
    GENERATED_BODY()
    FString ResidueName, ResidueSeq, Chain, RecordType;
};

USTRUCT()
struct FLigandMetadata {
    GENERATED_BODY()
    FString LigandName;
};

USTRUCT()
struct FResidueInfo
{
    GENERATED_BODY()
    FString RecordType, Chain, ResidueName, ResidueSeq;
    TArray<UStaticMeshComponent*> AtomMeshes, BondMeshes;
    // Store UNSCALED raw atom positions and element symbols
    TArray<FVector> AtomPositions;
    TArray<FString> AtomElements;
    TArray<TPair<int32, int32>> BondPairs;  // Bond connectivity
    TArray<int32> BondOrders;                // Bond orders
    bool bIsVisible = true;
};

USTRUCT()
struct FLigandInfo
{
    GENERATED_BODY()
    FString LigandName;
    TArray<UStaticMeshComponent*> AtomMeshes, BondMeshes;
    TArray<FString> AtomElements; // Element symbol per atom (aligned with AtomMeshes)
    TArray<FVector> AtomPositions; // Store UNSCALED positions for accurate calculations
    TArray<TPair<int32, int32>> BondPairs; // Store bond connectivity
    TArray<int32> BondOrders; // Bond order for each bond (aligned with BondPairs)
    bool bIsVisible = false;
};

// TreeView Node Object (must be UObject for TreeView)
UCLASS(BlueprintType)
class SPHERES_API UPDBTreeNode : public UObject
{
    GENERATED_BODY()
    
public:
    TMap<FString, FLigandInfo*> LigandMap;
    TMap<FString, FResidueInfo*> ResidueMap;
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    FString DisplayName;
    
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    FString NodeKey;
    
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    bool bIsChain = false;
    
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    bool bIsVisible = true;
    
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    FString ChainID;
    
    void Initialize(const FString& InDisplayName, const FString& InNodeKey, bool bInIsChain, const FString& InChainID = TEXT(""))
    {
        DisplayName = InDisplayName;
        NodeKey = InNodeKey;
        bIsChain = bInIsChain;
        bIsVisible = true;
        ChainID = InChainID;
    }
};

// ListView Node Object for SDF Molecules
UCLASS(BlueprintType)
class SPHERES_API UPDBMoleculeNode : public UObject
{
    GENERATED_BODY()
    
public:
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    FString MoleculeName;
    
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    FString MoleculeKey;
    
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    bool bIsVisible = true;
    
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    int32 AtomCount = 0;
    
    UPROPERTY(BlueprintReadOnly, Category = "PDB Viewer")
    int32 BondCount = 0;
    
    void Initialize(const FString& InName, const FString& InKey, bool bInIsVisible, int32 InAtomCount, int32 InBondCount)
    {
        MoleculeName = InName;
        MoleculeKey = InKey;
        bIsVisible = bInIsVisible;
        AtomCount = InAtomCount;
        BondCount = InBondCount;
    }
};

UCLASS()
class SPHERES_API APDBViewer : public AActor
{
    GENERATED_BODY()

public:
    APDBViewer();
    
    DECLARE_DYNAMIC_MULTICAST_DELEGATE(FOnResiduesLoaded);
    UPROPERTY(BlueprintAssignable, Category = "PDB Viewer") FOnResiduesLoaded OnResiduesLoaded;
    DECLARE_DYNAMIC_MULTICAST_DELEGATE(FOnLigandsLoaded);
    UPROPERTY(BlueprintAssignable, Category = "PDB Viewer") FOnLigandsLoaded OnLigandsLoaded;
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void SaveStructureToFile(const FString& FilePath);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void LoadStructureFromFile(const FString& FilePath);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void ClearCurrentStructure();
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void OpenSaveDialog();
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void OpenLoadDialog();
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void ToggleResidueVisibility(const FString& ResidueKey);
    
    // New TreeView functions
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") TArray<UPDBTreeNode*> GetChainNodes();
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") TArray<UPDBTreeNode*> GetResidueNodesForChain(const FString& ChainID);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void ToggleChainVisibility(const FString& ChainID);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void ToggleNodeVisibility(UPDBTreeNode* Node);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void PopulateTreeView(class UTreeView* TreeView);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") TArray<UObject*> GetChildrenForNode(UPDBTreeNode* Node);
    
    // ListView functions for SDF molecules
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") TArray<UPDBMoleculeNode*> GetMoleculeNodes();
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void PopulateMoleculeListView(class UListView* ListView);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void ToggleMoleculeVisibility(const FString& MoleculeKey);
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") void ToggleMoleculeNodeVisibility(UPDBMoleculeNode* Node);
    
    // Hydrogen generation functions
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void AddExplicitHydrogens();
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void RemoveExplicitHydrogens();
    
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer")
    void ToggleHydrogens();
    
    // Debug functions
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer|Debug")
    void DebugPrintLigandInfo();

    UFUNCTION(BlueprintCallable, Category = "PDB Viewer|Debug")
    int32 GetHydrogenCount() const;
    
    // Legacy functions (kept for backwards compatibility)
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") TArray<FString> GetResidueList() const;
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") TArray<FString> GetLigandList() const;
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") FString GetResidueDisplayName(const FString& ResidueKey) const;
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer") FString GetLigandDisplayName(const FString& LigandKey) const;

    // Get the currently visible ligand (if any)
    FLigandInfo* GetVisibleLigandInfo() const;
    TMap<FString, FLigandInfo*> LigandMap;
    TMap<FString, FResidueInfo*> ResidueMap;

    // Optional MD control widget class that will be spawned on BeginPlay if set
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "UI")
    TSubclassOf<class UMDControlWidget> MDControlWidgetClass;

    UPROPERTY()
    class UMDControlWidget* MDControlWidgetInstance;

    // Debug: highlight severe overlaps found by MMGBSA
    UFUNCTION(BlueprintCallable, Category = "PDB Viewer|Debug")
    void HighlightOverlapAtoms(const TArray<FMMOverlapInfo>& Overlaps);

    UFUNCTION(BlueprintCallable, Category = "PDB Viewer|Debug")
    void ClearOverlapMarkers();

protected:
    UFUNCTION()
    void OnLigandsLoadedHandler();

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "PDB Viewer")
    bool bAutoGenerateHydrogens = true;

    UPROPERTY() TArray<UStaticMeshComponent*> OverlapMarkers;

    virtual void BeginPlay() override;
    
    UPROPERTY() UStaticMesh* SphereMeshAsset;
    UPROPERTY() UStaticMesh* CylinderMeshAsset;
    UPROPERTY() UMaterial* SphereMaterialAsset;
    UPROPERTY() TArray<UStaticMeshComponent*> AllAtomMeshes;
    UPROPERTY() TArray<UStaticMeshComponent*> AllBondMeshes;
    
    FString CurrentPDBContent, CurrentStructureID;

    TSet<FString> ChainIDs; // Track all chains in the structure
    bool bHydrogensVisible = true;
    
    void FetchAndDisplayStructure(const FString& PDB_ID);
    void FetchFileAsync(const FString& URL, TFunction<void(bool, const FString&)> Callback);
    void ParsePDB(const FString& FileContent);
    void ParseMMCIF(const FString& FileContent);
    void ParseSDF(const FString& FileContent);
    void CreateResiduesFromAtomData(const TMap<FString, TMap<FString, FVector>>& ResidueAtoms, const TMap<FString, FResidueMetadata>& Metadata);
    void DrawProteinBondsAndConnectivity(const TMap<FString, FVector>& AtomPositions, FResidueInfo* ResInfo);
    void GenerateHydrogensForResidue(FResidueInfo* ResInfo);
    void FetchLigandBondsForHETATM(const FString& Key, const FString& Name, const TMap<FString, FVector>& Pos, FLigandInfo* LigInfo);
    void ParseLigandCIFForLigand(const FString& FileContent, const TMap<FString, FVector>& AtomPositions, FLigandInfo* LigInfo);
    FString NormalizeAtomID(const FString& In) const;
    int32 ParseBondOrder(const FString& OrderStr) const;
    void DrawSphere(float X, float Y, float Z, const FLinearColor& Color, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    void DrawBond(const FVector& Start, const FVector& End, int32 Order, const FString& Element1, const FString& Element2, USceneComponent* Parent, TArray<UStaticMeshComponent*>& OutArray);
    FLinearColor GetElementColor(const FString& Element) const;
    void ClearResidueMap();
    void ClearLigandMap();
    bool ShowFileDialog(bool bSave, FString& OutFilePath);
    
    // Hydrogen generation helpers
    int32 AddHydrogensToLigand(FLigandInfo* LigInfo);
};