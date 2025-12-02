#include "MySlateToolbar.h"
#include "Widgets/SBoxPanel.h"
#include "Widgets/Input/SButton.h"
#include "Widgets/Text/STextBlock.h"
#include "Engine/Engine.h" // For GEngine logs

void SMySlateToolbar::Construct(const FArguments& InArgs)
{
    // ... (SHorizontalBox creation code is correct) ...
    ChildSlot
    [
        SNew(SHorizontalBox)

        // --- Zoom Button ---
        + SHorizontalBox::Slot()
        .Padding(5.0f)
        .AutoWidth()
        [
            SNew(SButton)
            .OnClicked(this, &SMySlateToolbar::OnZoomClicked)
            .ContentPadding(FMargin(10, 5))
            [
                SNew(STextBlock)
                .Text(FText::FromString(TEXT("Zoom")))
            ]
        ]

        // --- Pan Button ---
        + SHorizontalBox::Slot()
        .Padding(5.0f)
        .AutoWidth()
        [
            SNew(SButton)
            .OnClicked(this, &SMySlateToolbar::OnPanClicked)
            .ContentPadding(FMargin(10, 5))
            [
                SNew(STextBlock)
                .Text(FText::FromString(TEXT("Pan")))
            ]
        ]

        // --- Rotate Button ---
        + SHorizontalBox::Slot()
        .Padding(5.0f)
        .AutoWidth()
        [
            SNew(SButton)
            .OnClicked(this, &SMySlateToolbar::OnRotateClicked)
            .ContentPadding(FMargin(10, 5))
            [
                SNew(STextBlock)
                .Text(FText::FromString(TEXT("Rotate")))
            ]
        ]
    ];
}

FReply SMySlateToolbar::OnZoomClicked() // Corrected: Use SMySlateToolbar::
{
    if (GEngine) GEngine->AddOnScreenDebugMessage(-1, 2.f, FColor::Yellow, TEXT("Zoom Clicked!"));
    return FReply::Handled();
}

// FIX: Change FReply FReply::OnPanClicked() to FReply SMySlateToolbar::OnPanClicked()
FReply SMySlateToolbar::OnPanClicked() 
{
    if (GEngine) GEngine->AddOnScreenDebugMessage(-1, 2.f, FColor::Yellow, TEXT("Pan Clicked!"));
    return FReply::Handled();
}

// FIX: Change FReply FReply::OnRotateClicked() to FReply SMySlateToolbar::OnRotateClicked()
FReply SMySlateToolbar::OnRotateClicked()
{
    if (GEngine) GEngine->AddOnScreenDebugMessage(-1, 2.f, FColor::Yellow, TEXT("Rotate Clicked!"));
    return FReply::Handled();
}
