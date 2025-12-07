#include "LigandVisibilityWidget.h"
#include "Ligands.h"
#include "Widgets/Text/STextBlock.h"
#include "Widgets/Input/SButton.h"
#include "Widgets/Layout/SScrollBox.h"
#include "Widgets/Layout/SBox.h"
#include "Widgets/Layout/SBorder.h"
#include "Widgets/SBoxPanel.h"  // SVerticalBox and SHorizontalBox
#include "Widgets/Layout/SSpacer.h"

void FLigandVisibilityWidget::Construct(const FArguments& InArgs)
{
    Ligands = InArgs._Ligands;

    // Horizontal box for Ligand buttons
    TSharedPtr<SHorizontalBox> ButtonBox = SNew(SHorizontalBox);

    if (Ligands.IsValid())
    {
        // ✅ Rename local variable to avoid shadowing
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
                .Text(FText::FromString(ButtonText))
                .ButtonColorAndOpacity(this, &FLigandVisibilityWidget::GetLigandButtonColor, i)
                .OnClicked(this, &FLigandVisibilityWidget::OnLigandButtonClicked, i)
            ];

            LigandButtons.Add(Button);
        }
    }

    // ✅ Wrap in a vertical box so the bar stays at the top
    ChildSlot
    [
        SNew(SVerticalBox)

        + SVerticalBox::Slot()
        .AutoHeight()
        .VAlign(VAlign_Top)
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

        + SVerticalBox::Slot()
        .FillHeight(1.0f)
        [
            SNew(SSpacer)
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

        if (LigandList.IsValidIndex(LigandIndex) && LigandList[LigandIndex].AtomSpheres.Num() > 0)
        {
            bool bCurrentVisibility = LigandList[LigandIndex].AtomSpheres[0]->IsVisible();

            // ✅ Use correct function name from ALigands
            Ligands->ToggleLigandVisibility(LigandIndex, !bCurrentVisibility);

            // ✅ Refresh button color immediately
            if (LigandButtons.IsValidIndex(LigandIndex) && LigandButtons[LigandIndex].IsValid())
            {
                LigandButtons[LigandIndex]->Invalidate(EInvalidateWidgetReason::Paint);
            }
        }
    }

    return FReply::Handled();
}
