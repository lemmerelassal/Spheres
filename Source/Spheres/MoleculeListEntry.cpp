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

    ApplyStyling();
}

void UMoleculeListEntry::ApplyStyling()
{
    // Style the molecule name
    if (TextMoleculeName)
    {
        TextMoleculeName->SetColorAndOpacity(FLinearColor(0.9f, 0.9f, 0.9f));
    }

    // Style the atom and bond counts
    if (TextAtomCount)
    {
        TextAtomCount->SetColorAndOpacity(FLinearColor(0.7f, 0.7f, 0.7f));
    }

    if (TextBondCount)
    {
        TextBondCount->SetColorAndOpacity(FLinearColor(0.7f, 0.7f, 0.7f));
    }

    // Style the toggle button
    if (ButtonToggleVisibility)
    {
        ButtonToggleVisibility->SetColorAndOpacity(FLinearColor(0.15f, 0.3f, 0.6f));
        
        // Style the button text
        UTextBlock* ButtonText = Cast<UTextBlock>(ButtonToggleVisibility->GetChildAt(0));
        if (ButtonText)
        {
            ButtonText->SetColorAndOpacity(FLinearColor::White);
        }
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
    UpdateButtonText();
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
        // Toggle visibility in the viewer
        Viewer->ToggleMoleculeNodeVisibility(CurrentMoleculeNode);

        // Update button text based on the NEW state after toggling
        UpdateButtonText();
    }
}

void UMoleculeListEntry::UpdateButtonText()
{
    if (!CurrentMoleculeNode || !ButtonToggleVisibility)
        return;

    UTextBlock* ButtonText = Cast<UTextBlock>(ButtonToggleVisibility->GetChildAt(0));
    if (ButtonText)
    {
        FString ToggleText = CurrentMoleculeNode->bIsVisible ? TEXT("Hide") : TEXT("Show");
        ButtonText->SetText(FText::FromString(ToggleText));

        // Change button color based on visibility state
        if (CurrentMoleculeNode->bIsVisible)
        {
            ButtonToggleVisibility->SetColorAndOpacity(FLinearColor(0.2f, 0.5f, 0.2f)); // Green when visible
        }
        else
        {
            ButtonToggleVisibility->SetColorAndOpacity(FLinearColor(0.5f, 0.2f, 0.2f)); // Red when hidden
        }
    }
}