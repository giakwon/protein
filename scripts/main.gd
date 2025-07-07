extends Node2D

var chain = ['H', 'P', 'P', 'H', 'P', 'H']
var res_scene = preload("res://scenes/residue.tscn")

func _ready():
	make_chain()

func make_chain():
	for i in range(len(chain)):
		var res = res_scene.instantiate()
		res.type = chain[i] 
		res.position = Vector2(100 + i * 80, 200)
		add_child(res)
