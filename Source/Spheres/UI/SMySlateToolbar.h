#pragma once

#include "CoreMinimal.h"
#include "Widgets/SCompoundWidget.h"

class SMySlateToolbar : public SCompoundWidget
{
public:
    SLATE_BEGIN_ARGS(SMySlateToolbar) {}
    SLATE_END_ARGS()

    void Construct(const FArguments& InArgs)
    {
        ChildSlot
        [
            // Example content
            SNew(STextBlock).Text(FText::FromString("Toolbar"))
        ];
    }
};
