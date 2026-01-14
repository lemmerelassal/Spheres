// InteractionControlWidget.h
#pragma once

#include "CoreMinimal.h"
#include "Blueprint/UserWidget.h"
#include "PDBViewer.h"
#include "InteractionControlWidget.generated.h"

// Forward declarations
class UTextBlock;
class UButton;
class UCheckBox;
class UScrollBox;
class UVerticalBox;

/**
 * Widget for controlling molecular interaction visualization
 */
UCLASS()
class SPHERES_API UInteractionControlWidget : public UUserWidget
{
    GENERATED_BODY()
    
public:
    UInteractionControlWidget(const FObjectInitializer& ObjectInitializer);
    
    // UI Components - bind these in the UMG designer
    UPROPERTY(meta = (BindWidget))
    UButton* CalculateButton;
    
    UPROPERTY(meta = (BindWidget))
    UCheckBox* HBondCheckBox;
    
    UPROPERTY(meta = (BindWidget))
    UCheckBox* SaltBridgeCheckBox;
    
    UPROPERTY(meta = (BindWidget))
    UCheckBox* PiStackCheckBox;
    
    UPROPERTY(meta = (BindWidget))
    UCheckBox* HydrophobicCheckBox;
    
    UPROPERTY(meta = (BindWidget))
    UCheckBox* ProteinProteinCheckBox;
    
    UPROPERTY(meta = (BindWidget))
    UCheckBox* ProteinLigandCheckBox;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* StatusText;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* HBondCountText;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* SaltBridgeCountText;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* PiStackCountText;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* HydrophobicCountText;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* TotalCountText;
    
    UPROPERTY(meta = (BindWidget))
    UScrollBox* InteractionListBox;
    
    // Reference to the PDB Viewer
    UPROPERTY(BlueprintReadWrite, Category = "Interactions")
    APDBViewer* PDBViewerRef;
    
    // Settings
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Interactions")
    bool bAutoCalculateOnLoad;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Interactions")
    bool bShowDetailedList;
    
    // Public API
    UFUNCTION(BlueprintCallable, Category = "Interactions")
    void SetPDBViewer(APDBViewer* Viewer);
    
    UFUNCTION(BlueprintCallable, Category = "Interactions")
    void TriggerCalculation();
    
    UFUNCTION(BlueprintCallable, Category = "Interactions")
    void OnInteractionsCalculated();
    
    UFUNCTION(BlueprintCallable, Category = "Interactions")
    void RefreshDisplay();
    
    UFUNCTION(BlueprintCallable, BlueprintPure, Category = "Interactions")
    int32 GetInteractionCountByType(EInteractionType Type) const;
    
protected:
    virtual void NativeConstruct() override;
    virtual void NativeDestruct() override;
    
private:
    // State
    bool bStructureReady = false;
    bool bIsCalculating = false;
    
    // Internal methods
    APDBViewer* FindPDBViewer();
    
    UFUNCTION()
    void OnStructureLoaded();
    
    UFUNCTION()
    void OnCalculateButtonClicked();
    
    UFUNCTION()
    void OnHBondCheckBoxChanged(bool bIsChecked);
    
    UFUNCTION()
    void OnSaltBridgeCheckBoxChanged(bool bIsChecked);
    
    UFUNCTION()
    void OnPiStackCheckBoxChanged(bool bIsChecked);
    
    UFUNCTION()
    void OnHydrophobicCheckBoxChanged(bool bIsChecked);
    
    void UpdateAllCounts();
    void UpdateStatusText(const FString& Status, const FLinearColor& Color);
    void UpdateInteractionList();
    void PopulateDetailedList();
    void SetUIEnabled(bool bEnabled);
};
