import random

deck_card = {}

def prepare_deck():
    card = 1
    while card <= 52:
        deck_card[card] = card
        card += 1

def pick_card():
    global picked_card
    c = random.randint(1, 52)
    while deck_card.get(c) == -1:
        c = random.randint(1, 52)
    picked_card = c
    deck_card[c] = -1


class check_card():
    def __init__(self):
        self.picked_card = picked_card
        self.num = picked_card % 13

    def information(self):
        if (self.picked_card - 1) // 13 == 0:
            self.suit = '黑桃'
        elif (self.picked_card - 1) // 13 == 1:
            self.suit = '红桃'
        elif (self.picked_card - 1) // 13 == 2:
            self.suit = '梅花'
        else:
            self.suit = '方块'

        if self.num == 1:
            self.card_name = 'A'
            self.value = 1
        elif self.num == 11:
            self.card_name = 'J'
            self.value = 10
        elif self.num == 12:
            self.card_name = 'Q'
            self.value = 10
        elif self.num == 0:
            self.card_name = 'K'
            self.value = 10
        else:
            self.card_name = str(self.num)
            self.value = self.num

    def check_Ace(self):
        global player_value
        if self.value == 1:
            while True:
                try:
                    answer = eval(input("你有一张" + self.suit + self.card_name + ",你想要它的值为1还是11呢？[1/11]"))
                    if answer == 1:
                        self.value = 1
                    elif answer == 11:
                        self.value = 11
                    else:
                        raise
                    player_value += self.value
                    break
                except:
                    print("请输入数字'1'或是'11'")
        else:
            player_value += self.value

    def dealer_check_Ace(self):
        global dealer_value
        if self.value == 1:
            if dealer_value <= 10:
                self.value = 11
            else:
                self.value = 1
            dealer_value += self.value
        else:
            dealer_value += self.value

    def show_card(self):
        global card
        card = self.suit + self.card_name


def dealer_steps(input):
    pick_card()
    input = check_card()
    input.information()
    input.dealer_check_Ace()
    input.show_card()


def player_steps(input):
    pick_card()
    input = check_card()
    input.information()
    input.check_Ace()
    input.show_card()


def gamble():
    try:
        global bet
        bet = eval(input('请为本次游戏下注：[10~500000](默认5000)'))
        if bet < 10 or bet > 500000:
            raise
    except:
        print('输入值错误，以默认值下注！')
        bet = 5000


def lose():
    global money
    money -= bet
    print('你失去了{}币，你目前的资产为{}甲鱼币'.format(bet, money))
    f = open('savedata.txt', 'w')
    f.write(str(money))
    f.close()


def win():
    global money
    money += bet
    print('你赢得了{}币，你目前的资产为{}甲鱼币'.format(bet, money))
    f = open('savedata.txt', 'w')
    f.write(str(money))
    f.close()


answer3 = 'Y'

while answer3 == 'Y' or answer3 == 'y':
    no_money = 0
    try:
        f = open('savedata.txt', 'r')
        money = eval(f.read())
        f.close()
        if money < 0:
            print('都没钱了还玩啥呀:(')
            no_money = 1
    except:
        money = 50000

    if no_money == 1:
        break

    for i in range(1):
        print('欢迎来到21点')
        print('你目前的资产为{}甲鱼币'.format(money))
        gamble()
        prepare_deck()
        dealer_value = 0
        player_value = 0
        dealercard1 = None
        dealer_steps(dealercard1)
        print('庄家的第1张牌是:' + card + ',庄家的第2张牌不展示，庄家现在点数为：' + str(dealer_value))
        dealercard2 = None
        dealer_steps(dealercard2)
        dealercard2 = card
        playercard1 = None
        player_steps(playercard1)
        print('你的第1张牌是：' + card)
        playercard2 = None
        player_steps(playercard2)
        print('你的第2张牌是：' + card + '，你现在的点数为：' + str(player_value))
        if player_value == 21:
            print('黑杰克！你赢得双倍赌注！！！')
            bet = bet * 2
            win()
        for i in range(3, 6):
            global break_level
            break_level = 0
            while True:
                try:
                    answer2 = input('是否继续要牌[Y/N]：')
                    if answer2 != 'Y' and answer2 != 'y' and answer2 != 'N' and answer2 != 'n':
                        raise
                    break
                except:
                    print("请输入正确的字符")

            if answer2 == 'Y' or answer2 == 'y':
                exec('palyer_card{} = {}'.format(i, None))
                player_steps("player_card{}".format(i))
                print('你的第{}张牌是：'.format(i) + card + '，你现在的点数为：' + str(player_value))
                if player_value >= 22:
                    print('爆炸！你输了！！')
                    lose()
                    break_level = 1
                    break
                elif i == 5:
                    print('你已有5张牌，你赢了')
                    win()
                    break_level = 1
                    break
            elif answer2 == 'N' or answer2 == 'n':
                break
            if break_level == 1:
                break
        if break_level == 1:
            break
        print("庄家的第2张牌是：" + dealercard2 + '，庄家现在点数为' + str(dealer_value))

        for i in range(3, 6):
            break_level = 0
            if dealer_value <= 16:
                exec('dealer_card{} = {}'.format(i, None))
                dealer_steps("dealer_card{}".format(i))
                print('庄家的第{}张牌是：'.format(i) + card + '，庄家现在的点数为：' + str(dealer_value))
                if dealer_value >= 22:
                    print('庄家爆炸！你赢了！！！')
                    win()
                    break_level = 1
                    break
                elif i == 5:
                    print('庄家已有5张牌，你输了！！')
                    lose()
                    break_level = 1
                    break
            if break_level == 1:
                break
        if break_level == 1:
            break
        print('庄家停止要牌')
        if dealer_value > player_value:
            print('庄家点数更大，你输了！！！')
            lose()
        elif dealer_value < player_value:
            print('你的点数更大，你赢了！！！')
            win()
        else:
            print('点数相同，平手！！！')
    if money < 0:
        print('破产咯，请充值甲鱼币（或者改存档也行）')
    print('***************************************************************')
    while True:
        try:
            answer3 = input('再来一局？[Y/N]：')
            if answer3 != 'Y' and answer3 != 'y' and answer3 != 'N' and answer3 != 'n':
                raise
            break
        except:
            print("请输入正确的字符")
