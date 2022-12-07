# PretextView tutorial

There are many ways to view a Hi-C contact map, but today we are going to use PretextView. To download the PretextView desktop application, click [here](https://github.com/wtsi-hpag/PretextView/releases), and pick a release that is suitable for your laptop. 

## How to curate your assemblies

### Step 1: Tweaking your settings

 To start the curation process, open the PretextView desktop application, click **load map**, and navigate to where you saved your pretext files. 
 
 When you load your Hi-C contact map, this is what it should look like, depending on which colour scheme you have chosen. For this tutorial we will be using "Blue-Orange Divergent".  
 
 Here, each square represents a scaffold (which after curation will hopefully all be of chromosome length). The red line in the centre shows where the strongest contact signals are between the DNA sequences. 

 Go ahead and look at the extensions for your contact maps. These overlays makes it easier to figure out where the unplaced and wrongly oriented scaffolds are supposed to go. 

 For this tutorial, turn on the **Gaps** extension, and turn the **Gamma Min** and **Gamma Mid** sliders down to zero, and the **Gamma Max** slider all the way up. This will make it easier to see where there are gaps in the assembly, and increase the contrast so the contact signals will be easier to interpret. 

 ### Step 2: Moving scaffolds in PretextView

For this step you´ll need a computer mouse. Before you do your first edit, try to orient yourself in the Hi-C contact map. Pressing **U** removes the main menu. Don´t worry, pressing **U** again brings it right back if you want to change any settings. 

Try to move around in the contact map. Scrolling the mouse wheel zooms you in and out. You´ll notice that the zoom is not affected by where your cursor is, so how do we zoom in on one particular scaffold? By dragging the map and placing the area you want to look at in the centre of the screen. To drag the map, click and hold the right mouse button, and move the mouse around. Move around the map for a bit, and look at the scaffolds. Do some look fragmented? Move to the far right bottom corner of the contact map. Are there many smaller scaffolds there?

When your comfortable with this movement, let´s try to bring up some of the other menus! Pressing **E** activates "Edit mode". When this mode is active, the cursor changes. Do you see that the further away from the red diagonal you move, the larger an area of the scaffold is marked by green? This green indicator shows which part of the scaffold you are picking up. Try to make an edit! Cut out a chunk of a scaffold and move it someplace else in the contact map. Don´t worry about whether the edit is correct, we´ll delete all the test edits before we start the proper curation process. 

Move to the far right of the contact map. Are there any smaller, unplaced scaffolds with clear red contact signals that you think would go well in any of the larger scaffolds? Hover over them, press the space bar, and move the cursor without clicking the left button. When you have oriented yourself to where you want to place the unplaced scaffold, click the left cursor. If the piece fits best on the end of one of the larger squares, press **S** while in Edit mode to toggle the Snap function. Do you notice that the unplaced scaffold "travels" differently across the contact map, in a skipping motion? This lets you "snap" the unplaced scaffold in place at the end of the larger scaffolds. 

Now that you know how to move around and make edits in PretextView, delete all your edits with **Q** while in Edit mode, and press **U** to check to see that all your test edits are gone. You are now ready to edit the assembly for real!

### Step 3: Editing your TPF

When editing a Hi-C contact map, you need to make *the same* edits to the TPF-file you generated with the Rapid curation suite. 

Bring up your TPF in your command line window with the command:

```
nano filename.tpf
```

### Step 4: Painting your scaffolds


### Step 5: Finishing your assembly