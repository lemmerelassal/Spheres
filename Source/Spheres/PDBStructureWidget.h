// PDBStructureWidget.h
#pragma once

#include "CoreMinimal.h"
#include "Blueprint/UserWidget.h"
#include "PDBStructureWidget.generated.h"

class UTreeView;
class APDBViewer;
class UPDBTreeNode;

UCLASS()
class SPHERES_API UPDBStructureWidget : public UUserWidget
{
    GENERATED_BODY()

public:
    virtual void NativeConstruct() override;
    virtual void NativeDestruct() override;

    // Bind this to your TreeView in the Designer
    UPROPERTY(BlueprintReadWrite, meta = (BindWidget))
    UTreeView* StructureTreeView;

protected:
    UPROPERTY()
    APDBViewer* PDBViewerRef;

    // Called when structure loads
    UFUNCTION()
    void OnStructureLoaded();

    // TreeView delegate for getting children - must match UE signature
    void OnGetItemChildren(UObject* Item, TArray<UObject*>& OutChildren);
};