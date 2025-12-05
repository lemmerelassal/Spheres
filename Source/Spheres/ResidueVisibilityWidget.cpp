#include "ResidueVisibilityWidget.h"
#include "BlueSphere.h"
#include "Widgets/Text/STextBlock.h"
#include "Widgets/Input/SButton.h"
#include "Widgets/Layout/SScrollBox.h"
#include "Widgets/Layout/SSpacer.h"

void FResidueVisibilityWidget::Construct(const FArguments& InArgs)
{
    BlueSphere = InArgs._BlueSphere;

    // Create a horizontal box to hold the buttons (instead of SVerticalBox)
    TSharedPtr<SHorizontalBox> ButtonBox = SNew(SHorizontalBox); // A horizontal box to hold buttons

    // Ensure we have the residues loaded from BlueSphere
    if (BlueSphere.IsValid())
    {
        // Get the AtomSpheres array using the getter
        TArray<UStaticMeshComponent*>& AtomSpheres = BlueSphere->GetAtomSpheres();

        // Iterate through all residues (or atoms depending on the structure)
        int32 Index = 0;
        for (auto& Residue : AtomSpheres)
        {
            const FString ButtonText = FString::Printf(TEXT("Residue %d"), Index + 1);

            // Add a button for each residue in the horizontal box
            ButtonBox->AddSlot()
                .AutoWidth()  // This ensures that each button is as wide as its content (no stretching)
                [
                    SNew(SButton)
                    .Text(FText::FromString(ButtonText))
                    .OnClicked(this, &FResidueVisibilityWidget::OnResidueButtonClicked, Index)
                ];

            ++Index;
        }
    }

    // Now that all buttons are added, assign the ButtonBox to the ChildSlot of the SScrollBox
    ChildSlot
    [
        SNew(SScrollBox) // Creates a scrollable container
        + SScrollBox::Slot()
        .Padding(2)
        [
            ButtonBox.ToSharedRef() // Add the horizontal box with buttons to the scroll box
        ]
    ];
}

FReply FResidueVisibilityWidget::OnResidueButtonClicked(int32 ResidueIndex)
{
    // Handle the visibility toggle for each residue
    if (BlueSphere.IsValid())
    {
        // Get the AtomSpheres array using the getter
        TArray<UStaticMeshComponent*>& AtomSpheres = BlueSphere->GetAtomSpheres();

        bool bCurrentVisibility = AtomSpheres[ResidueIndex]->IsVisible();
        AtomSpheres[ResidueIndex]->SetVisibility(!bCurrentVisibility);
    }

    return FReply::Handled();
}
