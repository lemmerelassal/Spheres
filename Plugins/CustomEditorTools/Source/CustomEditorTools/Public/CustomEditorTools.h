// Copyright Epic Games, Inc. All Rights Reserved.

#pragma once

#include "Modules/ModuleManager.h"
#include "Delegates/Delegate.h"

class FToolBarBuilder;
class FMenuBuilder;

class FCustomEditorToolsModule : public IModuleInterface
{
public:

	/** IModuleInterface implementation */
	virtual void StartupModule() override;
	virtual void ShutdownModule() override;
	
	/** This function will be bound to Command. */
	void PluginButtonClicked();
	
    // Declare the two new public handler functions here
	void Button2Clicked();
	void Button3Clicked();
	
private:

	void RegisterMenus();

private:
	TSharedPtr<class FUICommandList> PluginCommands;
    
    // Add this line to store the handle for the startup delegate
	FDelegateHandle StartupCallbackHandle; 
};
