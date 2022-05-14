# -*- coding: utf-8 -*-
"""
Created on Thu May 12 17:42:08 2022

@author: warbi
"""

import pgzrun
import pygame
import pgzero
import numpy as np
from pgzero.builtins import Actor
from random import randint
#from kivy import clock

WIDTH = 800
HEIGHT = 600

CENTER_X = WIDTH/2
CENTER_Y = HEIGHT/2

move_list = []
display_list = []




score = 0
current_move = 0
count = 4
dance_length  = 4

say_dance = False
show_countdown = True
moves_complete = False
game_over = False


dancer = Actor("dancer-start")
dancer.pos = CENTER_X + 5,CENTER_Y - 40

up = Actor("up")
up.pos = CENTER_X, CENTER_Y + 110
right = Actor("right")
right.pos = CENTER_X + 60, CENTER_Y + 170
down = Actor("down")
down.pos = CENTER_X, CENTER_Y + 230
left = Actor("left")
left.pos = CENTER_X - 60, CENTER_Y + 170

def draw():
    global game_over, score,say_dance, count, show_countdown
        
    if not game_over:
        screen.clear()
        screen.blit("stage",(0,0))
        dancer.draw()
        up.draw()
        down.draw()
        left.draw()
        right.draw()
        screen.draw.text("score: "+ str(score),color = "black", topleft = (10,10))
        screen.draw.text("Dance Length: " + str(dance_length), color="black", topright=(WIDTH - 10,10))
        if say_dance:
            screen.draw.text("Dance!", color="black",topleft=(CENTER_X - 65, 150), fontsize=60)
        if show_countdown:
            screen.draw.text(str(count), color="black",topleft=(CENTER_X - 8, 150), fontsize=60)
    else:
        screen.clear()
        screen.fill((255,255,255))
        screen.draw.text("Final score: " + str(score), color="black", topleft=(CENTER_X - 70, 60), fontsize=40)
        screen.draw.text("GAME OVER!", color="black",topleft=(CENTER_X - 130, 20), fontsize=60)
        display_high_scores()
    return

def display_high_scores():
    screen.draw.text('HIGH SCORES', (350, 150), color='black')
    y = 175
    position = 1
    for high_score in scores:
        screen.draw.text(str(position) + '.  ' + high_score, (350, y),
                         color='black')
        y += 25
        position += 1


scores = []
def update_high_scores():
    global score, scores
    #temporary score variable prevents value loss so it can be called on game over screen
    tempscore = score
    filename = r"high-scores.txt"
    scores = []
    with open(filename, 'r') as file:
        line = file.readline()
        high_scores = line.split()
        for high_score in high_scores:
            if(score > int(high_score)):
                scores.append(str(score) + ' ')
                score = int(high_score)
            else:
                scores.append(str(high_score) + ' ')
    with open(filename, 'w') as file:
        for high_score in scores:
            file.write(high_score)
    score = tempscore

def update():
    global game_over, current_move, moves_complete, dance_length
    if not game_over:
        if moves_complete:
            dance_length += 1
            generate_moves()

            moves_complete = False
            current_move = 0
            music.pause()
            sounds.happy_children.play()
            clock.schedule(music.unpause,3)
    else:
        music.stop()

def update_dancer(move):
    global game_over
    if not game_over:
        if move == 0:
            up.image = "up-lit"
            dancer.image = "dancer-up"
            clock.schedule(reset_dancer, 0.5)
        elif move == 1:
            right.image = "right-lit"
            dancer.image = "dancer-right"
            clock.schedule(reset_dancer, 0.5)
        elif move == 2:
            down.image = "down-lit"
            dancer.image = "dancer-down"
            clock.schedule(reset_dancer, 0.5)
        else:
            left.image = "left-lit"
            dancer.image = "dancer-left"
            clock.schedule(reset_dancer, 0.5)
        
        return
    
def reset_dancer():
    global game_over
    if not game_over:
        dancer.image = "dancer-start"
        up.image = "up"
        right.image = "right"
        down.image = "down"
        left.image = "left"
    return   

def on_key_up(key):
    global score, game_over, move_list, current_move
    if key == keys.UP:
        update_dancer(0)
        if move_list[current_move] == 0:
            score = score + 1
            next_move()
        else:
            game_over = True
            update_high_scores()
    elif key == keys.RIGHT:
        update_dancer(1)
        if move_list[current_move] == 1:
            score = score + 1
            next_move()
        else:
            game_over = True
            update_high_scores()
    elif key == keys.DOWN:
        update_dancer(2)
        if move_list[current_move] == 2:
            score = score + 1
            next_move()
        else:
            game_over = True
            update_high_scores()
    elif key == keys.LEFT:
        update_dancer(3)
        if move_list[current_move] == 3:
            score = score + 1
            next_move()
        else:
            game_over = True
            update_high_scores()
    return

def display_moves():
    global move_list, display_list, dance_length
    global say_dance, show_countdown, current_move
    if display_list:
        this_move = display_list[0]
        display_list = display_list[1:]
        if this_move == 0:
            update_dancer(0)
            clock.schedule(display_moves, 1)
        elif this_move == 1:
            update_dancer(1)
            clock.schedule(display_moves, 1)
        elif this_move == 2:
            update_dancer(2)
            clock.schedule(display_moves, 1)
        else:
            update_dancer(3)
            clock.schedule(display_moves, 1)
    else:
        say_dance = True
        show_countdown = False
    return

def countdown():
    global count, game_over, show_countdown
    if count > 1:
        count = count - 1
        clock.schedule(countdown, 1)
    else:
        show_countdown = False
        display_moves()
    return 

def generate_moves():
    global move_list, dance_length, count
    global show_countdown, say_dance
    count = 4
    move_list = []
    say_dance = False
    for move in range(0, dance_length):
        rand_move = randint(0, 3)
        move_list.append(rand_move)
        display_list.append(rand_move)
    show_countdown = True
    countdown()
    #this line for troublshooting
    #print(move_list)
    return

def next_move():
    global dance_length, current_move, moves_complete
    if current_move < dance_length - 1:
        current_move = current_move + 1
    else:
        moves_complete = True
    return

generate_moves()
music.play("turbo-strawberry")
music.queue("vanishing-horizon")
pgzrun.go()
