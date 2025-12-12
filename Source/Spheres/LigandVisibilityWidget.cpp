#include "LigandVisibilityWidget.h"
#include "Ligands.h"
#include "Widgets/Text/STextBlock.h"
#include "Widgets/Input/SButton.h"
#include "Widgets/Layout/SScrollBox.h"
#include "Widgets/Layout/SBox.h"
#include "Widgets/Layout/SBorder.h"
#include "Widgets/SBoxPanel.h"  // SVerticalBox and SHorizontalBox
#include "Widgets/Layout/SSpacer.h"
#include "Fonts/SlateFontInfo.h"
#include "Styling/CoreStyle.h"

void FLigandVisibilityWidget::Construct(const FArguments& InArgs)
{
    Ligands = InArgs._Ligands;

    // Horizontal box for Ligand buttons
    TSharedPtr<SHorizontalBox> ButtonBox = SNew(SHorizontalBox);

    if (Ligands.IsValid())
    {
        // Rename local variable to avoid shadowing
        const TArray<FLigandData>& LigandList = Ligands->GetLigands();

        for (int32 i = 0; i < LigandList.Num(); ++i)
        {
            const FString ButtonText = FString::Printf(TEXT("%s"), *LigandList[i].LigandName);

            TSharedPtr<SButton> Button;

            ButtonBox->AddSlot()
            .AutoWidth()
            .VAlign(VAlign_Center)
            .Padding(4, 2)
            [
                SAssignNew(Button, SButton)
                .ButtonColorAndOpacity(this, &FLigandVisibilityWidget::GetLigandButtonColor, i)
                .OnClicked(this, &FLigandVisibilityWidget::OnLigandButtonClicked, i)
                [
                    SNew(STextBlock)
                    .Text(FText::FromString(ButtonText))
                    .Font(FSlateFontInfo(FCoreStyle::GetDefaultFont(), 16, FName("Regular")))
                ]
            ];

            LigandButtons.Add(Button);
        }

        // Hide all ligands by default
        for (int32 i = 0; i < LigandList.Num(); ++i)
        {
            Ligands->ToggleLigandVisibility(i, false);
        }
    }

    // Wrap in a vertical box so the bar stays at the bottom
    ChildSlot
    [
        SNew(SVerticalBox)

        // Spacer to push scroll box down
        + SVerticalBox::Slot()
        .FillHeight(1.0f)
        [
            SNew(SSpacer)
        ]

        // Scroll box with buttons at the bottom
        + SVerticalBox::Slot()
        .AutoHeight()
        .VAlign(VAlign_Bottom)
        [
            SNew(SBorder)
            .Padding(5)
            .BorderBackgroundColor(FLinearColor(0.f, 0.f, 0.f, 0.4f))
            [
                SNew(SBox)
                .HeightOverride(60.f)
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
        ]
    ];
}

FSlateColor FLigandVisibilityWidget::GetLigandButtonColor(int32 LigandIndex) const
{
    if (Ligands.IsValid())
    {
        const TArray<FLigandData>& LigandList = Ligands->GetLigands();

        if (LigandList.IsValidIndex(LigandIndex) && LigandList[LigandIndex].AtomSpheres.Num() > 0)
        {
            bool bVisible = LigandList[LigandIndex].AtomSpheres[0]->IsVisible();
            return bVisible ? FLinearColor(0.0f, 0.2f, 1.0f, 1.0f)   // Blue when visible
                            : FLinearColor(0.8f, 0.0f, 0.0f, 1.0f);  // Red when hidden
        }
    }

    return FLinearColor(0.15f, 0.15f, 0.15f, 0.9f); // Default gray
}

FReply FLigandVisibilityWidget::OnLigandButtonClicked(int32 LigandIndex)
{
    if (Ligands.IsValid())
    {
        const TArray<FLigandData>& LigandList = Ligands->GetLigands();

        if (LigandList.IsValidIndex(LigandIndex))
        {
            bool bCurrentVisibility = false;
            if (LigandList[LigandIndex].AtomSpheres.Num() > 0)
            {
                bCurrentVisibility = LigandList[LigandIndex].AtomSpheres[0]->IsVisible();
            }

            if (bCurrentVisibility)
            {
                // If clicked ligand is visible, hide it (result: none visible)
                Ligands->ToggleLigandVisibility(LigandIndex, false);
            }
            else
            {
                // Hide all ligands first, then show selected one (ensures only one visible)
                for (int32 i = 0; i < LigandList.Num(); ++i)
                {
                    Ligands->ToggleLigandVisibility(i, false);
                }

                Ligands->ToggleLigandVisibility(LigandIndex, true);
            }

            // Refresh all button colors
            for (int32 i = 0; i < LigandButtons.Num(); ++i)
            {
                if (LigandButtons[i].IsValid())
                {
                    LigandButtons[i]->Invalidate(EInvalidateWidgetReason::Paint);
                }
            }
        }
    }

    return FReply::Handled();
}
