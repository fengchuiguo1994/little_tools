# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 19:26:19 2018

@author: huang
"""

from itertools import cycle
import random
import sys
import pygame
from pygame.locals import *

FPS = 30 #>25
SCREENWIDTH = 288
SCREENHEIGHT = 512
PIPEGAPSIZE = 100
BASEY = SCREENHEIGHT * 0.79
IMAGES,SOUNDS,HITMASKS = {},{},{}

PLAYERS_LIST = (
    #红色小鸟
    (
        'assets/sprites/redbird-upflap.png',      
        'assets/sprites/redbird-midflap.png',     
        'assets/sprites/redbird-downflap.png',      
    ),
    #蓝色小鸟
    (
        'assets/sprites/bluebird-upflap.png',        #分享前还是先分享自己的Python学习交流群：666468218
        'assets/sprites/bluebird-midflap.png',       #群内不定时分享干货，包括最新的python企业案例学习资料和零基础入门教程
        'assets/sprites/bluebird-downflap.png',      #欢迎初学和进阶中的小伙伴入群学习交流
    ),
    #黄色小鸟
    (
        'assets/sprites/yellowbird-upflap.png',
        'assets/sprites/yellowbird-midflap.png',
        'assets/sprites/yellowbird-downflap.png',
    ),
)
#调用背景
BACKGROUNDS_LIST = (
    'assets/sprites/background-day.png',
    'assets/sprites/background-night.png',
)
#调用管道
PIPES_LIST = (
    'assets/sprites/pipe-green.png',
    'assets/sprites/pipe-red.png',
)
try:
    xrange
except NameError:
    xrange = range
def main():
    global SCREEN, FPSCLOCK
    pygame.init()
    FPSCLOCK = pygame.time.Clock()
    SCREEN = pygame.display.set_mode((SCREENWIDTH, SCREENHEIGHT))
    pygame.display.set_caption('Flappy Bird')
    # 评分显示
    IMAGES['numbers'] = (
        pygame.image.load('assets/sprites/0.png').convert_alpha(),
        pygame.image.load('assets/sprites/1.png').convert_alpha(),
        pygame.image.load('assets/sprites/2.png').convert_alpha(),
        pygame.image.load('assets/sprites/3.png').convert_alpha(),#更多资料和源码请加QQ：666468218
        pygame.image.load('assets/sprites/4.png').convert_alpha(),#小编等待着大神和小白的到来，让我们一起学习交流
        pygame.image.load('assets/sprites/5.png').convert_alpha(),
        pygame.image.load('assets/sprites/6.png').convert_alpha(),
        pygame.image.load('assets/sprites/7.png').convert_alpha(),
        pygame.image.load('assets/sprites/8.png').convert_alpha(),
        pygame.image.load('assets/sprites/9.png').convert_alpha()
    )
    #游戏结束
    IMAGES['gameover'] = pygame.image.load('assets/sprites/gameover.png').convert_alpha()
    #欢迎界面
    IMAGES['message'] = pygame.image.load('assets/sprites/message.png').convert_alpha()
    #界面基础
    IMAGES['base'] = pygame.image.load('assets/sprites/base.png').convert_alpha()
    # 调用声音
    if 'win' in sys.platform:
        soundExt = '.wav'
    else:
        soundExt = '.ogg'
    SOUNDS['die']    = pygame.mixer.Sound('assets/audio/die' + soundExt)
    SOUNDS['hit']    = pygame.mixer.Sound('assets/audio/hit' + soundExt)
    SOUNDS['point']  = pygame.mixer.Sound('assets/audio/point' + soundExt)
    SOUNDS['swoosh'] = pygame.mixer.Sound('assets/audio/swoosh' + soundExt)
    SOUNDS['wing']   = pygame.mixer.Sound('assets/audio/wing' + soundExt)
    while True:
        #随机选择背景
        randBg = random.randint(0, len(BACKGROUNDS_LIST) - 1)
        IMAGES['background'] = pygame.image.load(BACKGROUNDS_LIST[randBg]).convert()

        #选择随机的小鸟
        randPlayer = random.randint(0, len(PLAYERS_LIST) - 1)
        IMAGES['player'] = (
            pygame.image.load(PLAYERS_LIST[randPlayer][0]).convert_alpha(),#分享前还是先分享自己的Python学习交流群：666468218
            pygame.image.load(PLAYERS_LIST[randPlayer][1]).convert_alpha(),#群内不定时分享干货，包括最新的python企业案例学习资料和零基础入门教程
            pygame.image.load(PLAYERS_LIST[randPlayer][2]).convert_alpha(),#欢迎初学和进阶中的小伙伴入群学习交流
        )
        #管道随机
        pipeindex = random.randint(0, len(PIPES_LIST) - 1)
        IMAGES['pipe'] = (
            pygame.transform.rotate(
                pygame.image.load(PIPES_LIST[pipeindex]).convert_alpha(), 180),
            pygame.image.load(PIPES_LIST[pipeindex]).convert_alpha(),
        )
        HITMASKS['pipe'] = (
            getHitmask(IMAGES['pipe'][0]),
            getHitmask(IMAGES['pipe'][1]),
        )
        HITMASKS['player'] = (
            getHitmask(IMAGES['player'][0]),
            getHitmask(IMAGES['player'][1]),
            getHitmask(IMAGES['player'][2]),
        )
        movementInfo = showWelcomeAnimation()
        crashInfo = mainGame(movementInfo)
        showGameOverScreen(crashInfo)
def showWelcomeAnimation():
    #显示欢迎界面动画的鸟
    playerIndex = 0
    playerIndexGen = cycle([0, 1, 2, 1])
    # 迭代器变化
    loopIter = 0
    playerx = int(SCREENWIDTH * 0.2)
    playery = int((SCREENHEIGHT - IMAGES['player'][0].get_height()) / 2)
    messagex = int((SCREENWIDTH - IMAGES['message'].get_width()) / 2)
    messagey = int(SCREENHEIGHT * 0.12)
    basex = 0
    baseShift = IMAGES['base'].get_width() - IMAGES['background'].get_width()
    #欢迎屏幕上的上下运动（额看到这里的，别想多了，我是一个正直的人）
    playerShmVals = {'val': 0, 'dir': 1}
    while True:
        for event in pygame.event.get():
            if event.type == QUIT or (event.type == KEYDOWN and event.key == K_ESCAPE):
                pygame.quit()
                sys.exit()
            if event.type == KEYDOWN and (event.key == K_SPACE or event.key == K_UP):
                SOUNDS['wing'].play()
                return {
                    'playery': playery + playerShmVals['val'],
                    'basex': basex,
                    'playerIndexGen': playerIndexGen,
                }
        if (loopIter + 1) % 5 == 0:
            playerIndex = next(playerIndexGen)
        loopIter = (loopIter + 1) % 30
        basex = -((-basex + 4) % baseShift)
        playerShm(playerShmVals)
        #绘制
        SCREEN.blit(IMAGES['background'], (0,0))
        SCREEN.blit(IMAGES['player'][playerIndex],
                    (playerx, playery + playerShmVals['val']))
        SCREEN.blit(IMAGES['message'], (messagex, messagey))
        SCREEN.blit(IMAGES['base'], (basex, BASEY))
        pygame.display.update()
        FPSCLOCK.tick(FPS)
def mainGame(movementInfo):
    score = playerIndex = loopIter = 0
    playerIndexGen = movementInfo['playerIndexGen']
    playerx, playery = int(SCREENWIDTH * 0.2), movementInfo['playery']
    basex = movementInfo['basex']
    baseShift = IMAGES['base'].get_width() - IMAGES['background'].get_width()
    #获得2个新的管道添加到上管，下管列表
    newPipe1 = getRandomPipe()
    newPipe2 = getRandomPipe()
    #上管列表
    upperPipes = [
        {'x': SCREENWIDTH + 200, 'y': newPipe1[0]['y']},
        {'x': SCREENWIDTH + 200 + (SCREENWIDTH / 2), 'y': newPipe2[0]['y']},
    ]
    #下管列表
    lowerPipes = [
        {'x': SCREENWIDTH + 200, 'y': newPipe1[1]['y']},
        {'x': SCREENWIDTH + 200 + (SCREENWIDTH / 2), 'y': newPipe2[1]['y']},
    ]
    pipeVelX = -4
    # 小鸟的速度，最大速度，向下的加速度，向上的加速度
    playerVelY    =  -9   # 小鸟的速度沿y轴，默认小鸟的拍打相同
    playerMaxVelY =  10   # 小鸟最大速度沿Y轴，最大下降速度
    playerMinVelY =  -8   # 小鸟最小速度沿y轴，最大提升速度
    playerAccY    =   1   # 小鸟向下的加速度
    playerRot     =  45   # 小鸟的轮换
    playerVelRot  =   3   # 角速度
    playerRotThr  =  20   
    playerFlapAcc =  -9   
    playerFlapped = False 
    while True:
        for event in pygame.event.get():
            if event.type == QUIT or (event.type == KEYDOWN and event.key == K_ESCAPE):
                pygame.quit()
                sys.exit()
            if event.type == KEYDOWN and (event.key == K_SPACE or event.key == K_UP):
                if playery > -2 * IMAGES['player'][0].get_height():
                    playerVelY = playerFlapAcc
                    playerFlapped = True
                    SOUNDS['wing'].play()
        # 坚持是否崩溃
        crashTest = checkCrash({'x': playerx, 'y': playery, 'index': playerIndex},
                               upperPipes, lowerPipes)
        if crashTest[0]:
            return {
                'y': playery,
                'groundCrash': crashTest[1],
                'basex': basex,
                'upperPipes': upperPipes,
                'lowerPipes': lowerPipes,
                'score': score,
                'playerVelY': playerVelY,
                'playerRot': playerRot
            }
        # 评分
        playerMidPos = playerx + IMAGES['player'][0].get_width() / 2
        for pipe in upperPipes:
            pipeMidPos = pipe['x'] + IMAGES['pipe'][0].get_width() / 2
            if pipeMidPos <= playerMidPos < pipeMidPos + 4:
                score += 1
                SOUNDS['point'].play()
        if (loopIter + 1) % 3 == 0:
            playerIndex = next(playerIndexGen)
        loopIter = (loopIter + 1) % 30
        basex = -((-basex + 100) % baseShift)
        # 旋转的小鸟
        if playerRot > -90:
            playerRot -= playerVelRot
        # 小鸟的运动
        if playerVelY < playerMaxVelY and not playerFlapped:
            playerVelY += playerAccY
        if playerFlapped:
            playerFlapped = False
            playerRot = 45
        playerHeight = IMAGES['player'][playerIndex].get_height()
        playery += min(playerVelY, BASEY - playery - playerHeight)
        # 向左边移动管道
        for uPipe, lPipe in zip(upperPipes, lowerPipes):
            uPipe['x'] += pipeVelX
            lPipe['x'] += pipeVelX
        # 当第一个管道即将接触屏幕左侧时添加新管道。
        if 0 < upperPipes[0]['x'] < 5:
            newPipe = getRandomPipe()
            upperPipes.append(newPipe[0])
            lowerPipes.append(newPipe[1])
        # 如果没有屏幕，移除第一个管道。
        if upperPipes[0]['x'] < -IMAGES['pipe'][0].get_width():
            upperPipes.pop(0)
            lowerPipes.pop(0)
        SCREEN.blit(IMAGES['background'], (0,0))
        for uPipe, lPipe in zip(upperPipes, lowerPipes):            #分享前还是先分享自己的Python学习交流群：666468218
            SCREEN.blit(IMAGES['pipe'][0], (uPipe['x'], uPipe['y']))#群内不定时分享干货，包括最新的python企业案例学习资料和零基础入门教程
            SCREEN.blit(IMAGES['pipe'][1], (lPipe['x'], lPipe['y']))#欢迎初学和进阶中的小伙伴入群学习交流
        SCREEN.blit(IMAGES['base'], (basex, BASEY))
        # 打印分数使玩家累积得分
        showScore(score)
        visibleRot = playerRotThr
        if playerRot <= playerRotThr:
            visibleRot = playerRot
        playerSurface = pygame.transform.rotate(IMAGES['player'][playerIndex], visibleRot)
        SCREEN.blit(playerSurface, (playerx, playery))
        pygame.display.update()
        FPSCLOCK.tick(FPS)
def showGameOverScreen(crashInfo):
    #死亡后的小鸟显示的游戏图像
    score = crashInfo['score']
    playerx = SCREENWIDTH * 0.2
    playery = crashInfo['y']
    playerHeight = IMAGES['player'][0].get_height()
    playerVelY = crashInfo['playerVelY']
    playerAccY = 2
    playerRot = crashInfo['playerRot']
    playerVelRot = 7
    basex = crashInfo['basex']
    upperPipes, lowerPipes = crashInfo['upperPipes'], crashInfo['lowerPipes']
    #撞击的声音和死亡的声音
    SOUNDS['hit'].play()
    if not crashInfo['groundCrash']:
        SOUNDS['die'].play()
    while True:
        for event in pygame.event.get():
            if event.type == QUIT or (event.type == KEYDOWN and event.key == K_ESCAPE):
                pygame.quit()
                sys.exit()
            if event.type == KEYDOWN and (event.key == K_SPACE or event.key == K_UP):
                if playery + playerHeight >= BASEY - 1:
                    return
        # 小鸟Y轴移动
        if playery + playerHeight < BASEY - 1:
            playery += min(playerVelY, BASEY - playery - playerHeight)
        # 小鸟的速度变化
        if playerVelY < 15:
            playerVelY += playerAccY
        if not crashInfo['groundCrash']:
            if playerRot > -90:
                playerRot -= playerVelRot
        SCREEN.blit(IMAGES['background'], (0,0))
        for uPipe, lPipe in zip(upperPipes, lowerPipes):
            SCREEN.blit(IMAGES['pipe'][0], (uPipe['x'], uPipe['y']))
            SCREEN.blit(IMAGES['pipe'][1], (lPipe['x'], lPipe['y']))
        SCREEN.blit(IMAGES['base'], (basex, BASEY))
        showScore(score)
        playerSurface = pygame.transform.rotate(IMAGES['player'][1], playerRot)
        SCREEN.blit(playerSurface, (playerx,playery))
        FPSCLOCK.tick(FPS)
        pygame.display.update()
def playerShm(playerShm):
    """oscillates the value of playerShm['val'] between 8 and -8"""
    if abs(playerShm['val']) == 8:
        playerShm['dir'] *= -1
    if playerShm['dir'] == 1:
         playerShm['val'] += 1
    else:
        playerShm['val'] -= 1
def getRandomPipe():
    #返回随机生成的管道
    # 上下管间隙
    gapY = random.randrange(0, int(BASEY * 0.6 - PIPEGAPSIZE))
    gapY += int(BASEY * 0.2)
    pipeHeight = IMAGES['pipe'][0].get_height()
    pipeX = SCREENWIDTH + 10
    return [
        {'x': pipeX, 'y': gapY - pipeHeight},  #上管
        {'x': pipeX, 'y': gapY + PIPEGAPSIZE}, #下管
    ]
def showScore(score):
    """displays score in center of screen"""
    scoreDigits = [int(x) for x in list(str(score))]
    totalWidth = 0 # 所有要打印的数字的总宽度
    for digit in scoreDigits:
        totalWidth += IMAGES['numbers'][digit].get_width()
    Xoffset = (SCREENWIDTH - totalWidth) / 2
    for digit in scoreDigits:
        SCREEN.blit(IMAGES['numbers'][digit], (Xoffset, SCREENHEIGHT * 0.1))
        Xoffset += IMAGES['numbers'][digit].get_width()
def checkCrash(player, upperPipes, lowerPipes):
    pi = player['index']
    player['w'] = IMAGES['player'][0].get_width()
    player['h'] = IMAGES['player'][0].get_height()
    # 如果玩家撞到地面
    if player['y'] + player['h'] >= BASEY - 1:
        return [True, True]
    else:
        playerRect = pygame.Rect(player['x'], player['y'],
                      player['w'], player['h'])
        pipeW = IMAGES['pipe'][0].get_width()
        pipeH = IMAGES['pipe'][0].get_height()
        for uPipe, lPipe in zip(upperPipes, lowerPipes):
            # 上、下管矩形
            uPipeRect = pygame.Rect(uPipe['x'], uPipe['y'], pipeW, pipeH)#分享前还是先分享自己的Python学习交流群：666468218
            lPipeRect = pygame.Rect(lPipe['x'], lPipe['y'], pipeW, pipeH)#群内不定时分享干货，包括最新的python企业案例学习资料和零基础入门教程
            pHitMask = HITMASKS['player'][pi]                            #欢迎初学和进阶中的小伙伴入群学习交流
            uHitmask = HITMASKS['pipe'][0]
            lHitmask = HITMASKS['pipe'][1]
            # 如果鸟与管道相撞
            uCollide = pixelCollision(playerRect, uPipeRect, pHitMask, uHitmask)
            lCollide = pixelCollision(playerRect, lPipeRect, pHitMask, lHitmask)
            if uCollide or lCollide:
                return [True, False]
    return [False, False]
def pixelCollision(rect1, rect2, hitmask1, hitmask2):
    #检查两个对象是否碰撞
    rect = rect1.clip(rect2)
    if rect.width == 0 or rect.height == 0:
        return False
    x1, y1 = rect.x - rect1.x, rect.y - rect1.y
    x2, y2 = rect.x - rect2.x, rect.y - rect2.y
    for x in xrange(rect.width):
        for y in xrange(rect.height):
            if hitmask1[x1+x][y1+y] and hitmask2[x2+x][y2+y]:
                return True
    return False
def getHitmask(image):
    mask = []
    for x in xrange(image.get_width()):
        mask.append([])
        for y in xrange(image.get_height()):
            mask[x].append(bool(image.get_at((x,y))[3]))
    return mask
if __name__ == '__main__':
    main()
