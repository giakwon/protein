class_name Residue
extends Area2D

@export var type: String = "H"
var is_dragging = false

#@onready var sprite = $Circle
@onready var label = $Label

func _ready():
	label.text = type
	if type == "H":
		$ColorRect.color = Color.DARK_ORANGE
	else:
		$ColorRect.color = Color.BLUE
		
#
#func _input_event(viewport, event, shape_idx):
	#if event is InputEventMouseButton and event.button_index == MOUSE_BUTTON_LEFT:
		#is_dragging = event.pressed
#
#func _process(_delta):
	#if is_dragging:
		#global_position = get_global_mouse_position()
	#elif not Input.is_mouse_button_pressed(MOUSE_BUTTON_LEFT) and is_dragging:
		#is_dragging = false
		#global_position = snap_to_grid(global_position)
#
#func snap_to_grid(pos: Vector2) -> Vector2:
	#var grid_size = 64
	#return Vector2(round(pos.x / grid_size), round(pos.y / grid_size)) * grid_size
