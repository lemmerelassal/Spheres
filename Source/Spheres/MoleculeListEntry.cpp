// MoleculeListEntry.cpp
#include "MoleculeListEntry.h"
#include "PDBViewer.h"
#include "Components/TextBlock.h"
#include "Components/Button.h"
#include "Kismet/GameplayStatics.h"

void UMoleculeListEntry::NativeConstruct()
{
    Super::NativeConstruct();

    if (ButtonToggleVisibility)
    {
        ButtonToggleVisibility->OnClicked.AddDynamic(this, &UMoleculeListEntry::OnToggleClicked);
    }
}

void UMoleculeListEntry::NativeOnListItemObjectSet(UObject* ListItemObject)
{
    // This is called when the ListView assigns a molecule node to this entry
    CurrentMoleculeNode = Cast<UPDBMoleculeNode>(ListItemObject);

    if (!CurrentMoleculeNode)
        return;

    // Update the UI elements
    if (TextMoleculeName)
    {
        TextMoleculeName->SetText(FText::FromString(CurrentMoleculeNode->MoleculeName));
    }

    if (TextAtomCount)
    {
        FString AtomText = FString::Printf(TEXT("Atoms: %d"), CurrentMoleculeNode->AtomCount);
        TextAtomCount->SetText(FText::FromString(AtomText));
    }

    if (TextBondCount)
    {
        FString BondText = FString::Printf(TEXT("Bonds: %d"), CurrentMoleculeNode->BondCount);
        TextBondCount->SetText(FText::FromString(BondText));
    }

    // Update button appearance based on visibility
    if (ButtonToggleVisibility)
    {
        UTextBlock* ButtonText = Cast<UTextBlock>(ButtonToggleVisibility->GetChildAt(0));
        if (ButtonText)
        {
            FString ToggleText = CurrentMoleculeNode->bIsVisible ? TEXT("Hide") : TEXT("Show");
            ButtonText->SetText(FText::FromString(ToggleText));
        }
    }
}

void UMoleculeListEntry::OnToggleClicked()
{
    if (!CurrentMoleculeNode)
        return;

    // Find the PDBViewer actor in the world
    APDBViewer* Viewer = Cast<APDBViewer>(
        UGameplayStatics::GetActorOfClass(GetWorld(), APDBViewer::StaticClass())
    );

    if (Viewer)
    {
        Viewer->ToggleMoleculeNodeVisibility(CurrentMoleculeNode);

        // Update button text
        if (ButtonToggleVisibility)
        {
            UTextBlock* ButtonText = Cast<UTextBlock>(ButtonToggleVisibility->GetChildAt(0));
            if (ButtonText)
            {
                FString ToggleText = CurrentMoleculeNode->bIsVisible ? TEXT("Hide") : TEXT("Show");
                ButtonText->SetText(FText::FromString(ToggleText));
            }
        }
    }
}
