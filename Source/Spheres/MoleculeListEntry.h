// MoleculeListEntry.h
#pragma once

#include "CoreMinimal.h"
#include "Blueprint/UserWidget.h"
#include "Blueprint/IUserObjectListEntry.h"
#include "MoleculeListEntry.generated.h"

class UTextBlock;
class UButton;
class UPDBMoleculeNode;

UCLASS()
class SPHERES_API UMoleculeListEntry : public UUserWidget, public IUserObjectListEntry
{
    GENERATED_BODY()

public:
    // IUserObjectListEntry interface
    virtual void NativeOnListItemObjectSet(UObject* ListItemObject) override;

protected:
    virtual void NativeConstruct() override;

    // Bind these in the UMG Designer (BindWidget)
    UPROPERTY(meta = (BindWidget))
    UTextBlock* TextMoleculeName;

    UPROPERTY(meta = (BindWidget))
    UTextBlock* TextAtomCount;

    UPROPERTY(meta = (BindWidget))
    UTextBlock* TextBondCount;

    UPROPERTY(meta = (BindWidget))
    UButton* ButtonToggleVisibility;

private:
    UFUNCTION()
    void OnToggleClicked();

    UPROPERTY()
    UPDBMoleculeNode* CurrentMoleculeNode;
};
