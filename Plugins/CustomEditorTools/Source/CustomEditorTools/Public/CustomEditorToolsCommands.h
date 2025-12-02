// Copyright Epic Games, Inc. All Rights Reserved.

#pragma once

#include "Framework/Commands/Commands.h"
#include "CustomEditorToolsStyle.h"

class FCustomEditorToolsCommands : public TCommands<FCustomEditorToolsCommands>
{
public:
	FCustomEditorToolsCommands()
		: TCommands<FCustomEditorToolsCommands>(TEXT("CustomEditorTools"), NSLOCTEXT("Contexts", "CustomEditorTools", "CustomEditorTools Plugin"), NAME_None, FCustomEditorToolsStyle::GetStyleSetName())
	{
	}

	// TCommands<> interface
	virtual void RegisterCommands() override;

public:
	TSharedPtr<FUICommandInfo> PluginAction;

	// Declare the new commands here
	TSharedPtr<FUICommandInfo> Button2Action;
	TSharedPtr<FUICommandInfo> Button3Action;
};
