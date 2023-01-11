/**
* Name: button
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model button

/**
 * A generic button interface that allows to select items and to move selected items.
 */
global {
	/**
	 * The button currently under mouse.
	 */
	Button current_button_focus;
	/**
	 * The last button that catched a click event.
	 */
	Button last_clicked_button;
	/**
	 * The item eventually selected by a clicked button.
	 */
	agent selected_item;
	
	/**
	 * Used to disable the mouse_up events that immediatly follows the first mouse_down event on a button.
	 */
	bool disable_up <- false;
	
	/**
	 * Find buttons among the specified handled_buttons species that are under the
	 * cursor, and handles the move of currently selected item under the mouse cursor.
	 * 
	 * The mouse_enter and mouse_leave actions are called on the button under the
	 * cursor when the mouse enters or leaves the button area, only if there is
	 * currently no selected item.
	 */
	action mouse_move_buttons_list(list<species<Button>> handled_buttons) {
		Button button_under_mouse <- nil;
		
		int i <- 0;
		loop while: ((button_under_mouse = nil) and (i < length(handled_buttons))) {
			list<Button> buttons_under_mouse <- handled_buttons[i] overlapping #user_location;
			if length(buttons_under_mouse) > 0 {
				button_under_mouse <- buttons_under_mouse[0];
			}
			i <- i+1;
		}
			
		// Mouse leave is triggered if a button was previously in focus and no button is
		// currently in focus or the button currently in focus is not the same as the previous one.
		if current_button_focus != nil and (button_under_mouse = nil or button_under_mouse != current_button_focus) {
			ask current_button_focus {
				do mouse_leave;
			}
			current_button_focus <- nil;
		}	

		// Mouse enter is not triggered if an item is currently selected
		if (selected_item = nil) {
			// The mouse enter action is called only when the mouse enters a new button
			// (the action is not triggered when the mouse moves within a button)
			if button_under_mouse != nil and button_under_mouse != current_button_focus {
				current_button_focus <- button_under_mouse;
				ask current_button_focus {
					do mouse_enter;
				}
			}
		} else {
			// An item is selected, so move the item but do not enters any other button
			// until the item is released.
			selected_item.location <- #user_location;
		}

	}
	
	/**
	 * Action triggered when a button is clicked and no item is currently selected.
	 */
	action mouse_down_buttons {
		if current_button_focus != nil and selected_item = nil {
			ask current_button_focus {
				last_clicked_button <- self;
				selected_item <- self.click();
				// Disables the mouse_up event just following this mouse_down
				disable_up <- true;
			}
		}
	}
	
		
	action mouse_up_buttons {
		if(disable_up) {
			// Does nothing this time, but handles the next mouse_up
			disable_up <- false;
		} else if (last_clicked_button != nil) {
			ask last_clicked_button {
				do post_click;
			}
			last_clicked_button <- nil;
			selected_item <- nil;
		}
	}
}

/**
 * A Button interface. Mouse events are catched within the button shape, that can be overidden has needed.
 */
species Button {
	/**
	 * Base button size. The default shape of the button is a circle of diameter button_size.
	 */
	float button_size;
	
	init {
		shape <- circle(0.5*button_size); // 0.5 factor for radius
	}
	
	/**
	 * Action called when the mouse enters the button.
	 */
	action mouse_enter virtual: true;
	/**
	 * Action called when the mouse leaves the button.
	 */
	action mouse_leave virtual: true;
	/**
	 * Action called when the button is clicked.
	 * The action is called on the mouse_down event.
	 * 
	 * If the action returns an agent, the agent becomes the selected_item and moves
	 * with the cursor.
	 * 
	 * If no item selection is needed, the action can safely return nil.
	 */
	action click virtual: true type: agent;
	/**
	 * Action called on the next click after the button was previously clicked.
	 * The action is called even if the second click does not occur on the button.
	 * This can be useful to release/delete/move an item selected by the click action.
	 * 
	 * The action is called on the mouse_up event, so that the selected item can be
	 * safely used on the mouse_down event of the second click and die in the post_click
	 * action.
	 */
	action post_click virtual: true;
	
}
