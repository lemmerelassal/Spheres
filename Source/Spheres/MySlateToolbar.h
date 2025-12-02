#pragma once

#include "Widgets/SCompoundWidget.h"

class SMySlateToolbar : public SCompoundWidget
{
public:
    SLATE_BEGIN_ARGS(SMySlateToolbar) {}
    SLATE_END_ARGS()

    void Construct(const FArguments& InArgs);

protected:
    // Handlers for button clicks
    FReply OnZoomClicked();
    FReply OnPanClicked();
    FReply OnRotateClicked();
};
