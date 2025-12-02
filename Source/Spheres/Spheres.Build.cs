// Copyright Epic Games, Inc. All Rights Reserved.

using UnrealBuildTool;

public class Spheres : ModuleRules
{
	public Spheres(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;

        if (Target.Type == TargetType.Editor)
        {
            PrivateDependencyModuleNames.AddRange(new string[] { "LevelEditor", "Slate", "SlateCore", "EditorStyle" });
        }
	
		PublicDependencyModuleNames.AddRange(new string[] { "Core", "CoreUObject", "Engine", "InputCore",  "EnhancedInput",
            "Json",
            "JsonUtilities",
			
            "Slate",
            "SlateCore",
			"UMG" });

		PrivateDependencyModuleNames.AddRange(new string[] { "Eigen" });

		// Uncomment if you are using Slate UI
		// PrivateDependencyModuleNames.AddRange(new string[] { "Slate", "SlateCore" });
		
		// Uncomment if you are using online features
		// PrivateDependencyModuleNames.Add("OnlineSubsystem");

		// To include OnlineSubsystemSteam, add it to the plugins section in your uproject file with the Enabled attribute set to true
	}
}
