// =========================================================================
// MoleculeListWidget.h
#pragma once

#include "CoreMinimal.h"
#include "Blueprint/UserWidget.h"
#include "MoleculeListWidget.generated.h"

class UListView;
class APDBViewer;

UCLASS()
class SPHERES_API UMoleculeListWidget : public UUserWidget
{
    GENERATED_BODY()

public:
    // Call this to populate the list with molecules from a PDBViewer
    UFUNCTION(BlueprintCallable, Category = "Molecule List")
    void PopulateList(APDBViewer* Viewer);

    // Automatically find and populate from first PDBViewer in world
    UFUNCTION(BlueprintCallable, Category = "Molecule List")
    void AutoPopulate();

protected:
    virtual void NativeConstruct() override;

    // Bind this in the UMG Designer (BindWidget)
    UPROPERTY(meta = (BindWidget))
    UListView* MoleculeListView;

private:
    UFUNCTION()
    void OnLigandsLoaded();

    UPROPERTY()
    APDBViewer* CachedViewer;
};

