#include "ResidueVisibilityWidget.h"
#include "BlueSphere.h"
#include "Widgets/Text/STextBlock.h"
#include "Widgets/Input/SButton.h"
#include "Widgets/Layout/SScrollBox.h"
#include "Widgets/Layout/SBox.h"
#include "Widgets/Layout/SBorder.h"

void FResidueVisibilityWidget::Construct(const FArguments& InArgs)
{
    BlueSphere = InArgs._BlueSphere;

    // A horizontal box that holds the residue buttons
    TSharedPtr<SHorizontalBox> ButtonBox = SNew(SHorizontalBox);

    if (BlueSphere.IsValid())
    {
        const TArray<FResidueData>& Residues = BlueSphere->GetResidues();

        for (int32 i = 0; i < Residues.Num(); ++i)
        {
            const FString ButtonText = FString::Printf(TEXT("%s"), *Residues[i].ResidueName);

            ButtonBox->AddSlot()
            .AutoWidth()                // Keep each button sized to content
            .VAlign(VAlign_Center)      // Vertically center buttons
            .Padding(4, 2)
            [
                SNew(SButton)
                .Text(FText::FromString(ButtonText))
                .ButtonColorAndOpacity(FLinearColor(0.15f, 0.15f, 0.15f, 0.9f))
                .OnClicked(this, &FResidueVisibilityWidget::OnResidueButtonClicked, i)
            ];
        }
    }

    // Wrap in a scroll box that only takes a small vertical slice at the top/bottom of the screen
    ChildSlot
    [
        SNew(SBorder)
        .Padding(5)
        .BorderBackgroundColor(FLinearColor(0.f, 0.f, 0.f, 0.4f))
        [
            SNew(SBox)
            .HeightOverride(60.f) // 👈 Limit the bar height so it doesn't fill the screen
            [
                SNew(SScrollBox)
                .Orientation(Orient_Horizontal)
                .ScrollBarVisibility(EVisibility::Collapsed)
                + SScrollBox::Slot()
                .Padding(2)
                [
                    ButtonBox.ToSharedRef()
                ]
            ]
        ]
    ];
}

FReply FResidueVisibilityWidget::OnResidueButtonClicked(int32 ResidueIndex)
{
    if (BlueSphere.IsValid())
    {
        const TArray<FResidueData>& Residues = BlueSphere->GetResidues();
        if (Residues.IsValidIndex(ResidueIndex) && Residues[ResidueIndex].AtomSpheres.Num() > 0)
        {
            bool bCurrentVisibility = Residues[ResidueIndex].AtomSpheres[0]->IsVisible();
            BlueSphere->ToggleResidueVisibility(ResidueIndex, !bCurrentVisibility);
        }
    }

    return FReply::Handled();
}
