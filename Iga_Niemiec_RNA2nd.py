# zadanie:
# zamiana struktury 2nd RNA zapisane w formacie dot-bracket na graf przedstawiający rozmieszczenie poszczególnych struktur

import re
in_file = open("C:/Users/Igusia/Documents/example_complicated.txt")
in_file2 = in_file.read().replace("\n", "")
input_file = str(in_file2)
#print("Plik wejsciowy: \n", input_file)
#print()
# importowanie pliku wejściowego, usunięcie newline i zamiana na string

bases = ["A", "C", "U", "G", "N"]
signs = [".", "(", ")", "[", "]", ">", "<", "{", "}"]
# listy zasad azotowych (po nich następuje właściwa sekwencja dot-bracket), oraz znaków w poszukiewanym formacie

def find_seq(input_seq):
    sequence = ""
    for i in range(len(input_seq)):
        if input_seq[i] in bases and input_seq[i+1] in signs and input_seq[i+2] in signs and input_seq[i+3] in signs:
            sequence = (input_seq[i+1:])
    return sequence
# funkcja zwracająca sekwencję w formacie dot-bracket z wczytywanego pliku
# szukająca pierwszych 3 znaków formatu d-b następujących po sekwencji zasad

dot_bracket_seq = find_seq(input_file)
#print("Sekwencja w formacie dot-bracket: \n", dot_bracket_seq)
# wypisanie sekwencji w formacie d-b

##########################################################################


def clean_sequence(seq):  # usuwa pseudowęzły oraz jednoniciowe fragmenty końcowe z sekwencji
    rev_seq = seq[::-1]
    left_onestrand = 0
    right_onestrand = 0
    for i in seq:
        if i != "." and i != "(" and i != ")":
            seq = seq.replace(i, ".")
            # usuwanie pseudowęzłów - zamiana wszystkich dodatkowych znaków na "."
    for i in range(len(seq)):
        if seq[i] == "(":
            seq = seq[i:]
            left_onestrand = i
            break
            # odnajdywanie 1szego "(" od lewej strony i odcinanie fragmentu po lewej stronie
    for j in range(len(rev_seq)):
        if rev_seq[j] == ")":
            seq = seq[:(len(seq)-j)] # len(seq) bo po odcięciu w poprzedniej pętli jej dł jest zmieniona
            right_onestrand = j
            break
            # analogicznie - po prawej stronie
    return seq, left_onestrand, right_onestrand

input_sequence = clean_sequence(dot_bracket_seq)[0]
left_onestrand = clean_sequence(dot_bracket_seq)[1]
right_onestrand = clean_sequence(dot_bracket_seq)[2]

print("Sekwencja wejściowa w formacie dot-bracket:", input_sequence)
print()

#################################################################################


def find_hairpins(seq): # znajduje struktury spinek do włosów w sekwencji, zwraca listy  ze słownikami
    hairpin_pat = "\(\.*\)" # wyrazenie regularne znajdujace spinki (dowolna ilość "." zakończona "(" z lewej i ")" z prawej strony)
    hairpin_pattern = re.compile(hairpin_pat) # kompilacja wyrażenia regularnego
    structures = []

    for m in hairpin_pattern.finditer(seq): # pętla dla każdego patternu wyszukanego w sekwencji
        structures.append([{'type': 'hairpin', 'start': m.start() + 1, 'end': m.end() - 1, 'length': (m.end()-1) - (m.start()+1)}])
        #  każda str będzie zapisana w osobnej liście skł się ze słowników z informacjami - na razie są tam inf o znalezionej spince
    return structures

structures = find_hairpins(input_sequence)
print("Znalezione struktury spinek do włosów:", structures)
print()

################################################################################

# zapisywanie sekwencji na prawo i na lewo od każdej spinki

def left_and_right_structure_side(seq, structures_list): # wypisuje części sekwencji wejściowej po prawej i lewej stronie od każdej struktury
    left_str_side = []
    right_str_side = []

    for n in range(len(structures_list)):
        left_side = seq[: structures_list[n][0]['start']]
        left_str_side.append(left_side) # lista sekwencji po lewej stronie każdej ze znalezionch struktur
        right_side = seq[structures_list[n][0]['end']:]
        right_str_side.append(right_side) # analogiczna lista z sekwencjami na prawo od spinek
    return left_str_side, right_str_side
    # bez uwzględnienia znaków wchodzących w skłąd spinki
    # stworzenie list sekwencji "prawostronnych" oraz "lewostronnnych"

left_structures_side = left_and_right_structure_side(input_sequence, structures)[0]
right_structures_side = left_and_right_structure_side(input_sequence, structures)[1]

print("Lista sekwencji po lewej stronie struktur:", left_structures_side)
print("Lista sekwencji po prawej stronie struktur:", right_structures_side)
print()

#################################################################################################
# funkcja odnajdująca pnie/pętle/bulge

def find_slb(left, right, struct):  # odnajduje pnie/pętle/bulge wew(stems/bulges/loops) połączone liniowo ze spinkami lub str wyższego rzędu
    # i dopisuje je do listy w postaci słowników z inf o strukturze

    ending_left_pat = "\.*\)" # wyrażenie regularne odnajdujące koniec lewostronny danej struktury liniowej (tj kropki i nawias w odwrotą stronę)
    ending_left_pattern = re.compile(ending_left_pat) # kompilacja wyrażenia regularnego
    ending_right_pat = "\.*\(" # jw. - koniec prawostronny
    ending_right_pattern = re.compile(ending_right_pat)


    for n in range(len(left)):  # zakres odp ilości znalezionych struktur
        rev_left = left[n][::-1] # odwrócenie sekwencji "lewych cześci spinek" aby łatwiejsze było porównywanie sekwencji "prawej" i "lewej"
        sequence_in_progress = 0  # zmienna oznaczająca strukturę (pień/pętlę/bulge) w której pętli obecnie znajduje się funkcja
        i = 0 # zmienna oznaczająca znak w sekwencji "lewej"
        j = 0 # zmienna oznaczająca znak w sekwencji "prawej"

        ending_places_left = [] # lista oznaczająca znalezione na podstawie patternu miejsca kończące funkcji po lewej stronie
        ending_places_right = [] # j.w - po prawej stronie

        for end_l in ending_left_pattern.finditer(rev_left): # pętla dla każdego patternu "kończącego" wyszukanego w sekwencji "lewostronnej"
            ending_places_left.append(end_l.start()) # dodanie indeksu znaku rozpoczynającego pattern do listy
        for end_r in ending_right_pattern.finditer(right[n]): # j.w - dotyczy strony prawej
            ending_places_right.append(end_r.start())

        if len(ending_places_left) == 0: # jeśli nie znaleziono żadnego patternu "konczącego" funkcja będzie działać aż do końca sekwencji
            end_place_left = len(rev_left)
        else:
            end_place_left = ending_places_left[0] # jeśli zostanie znaleziony pattern "kończący" funkcja będzie działać aż do 1szego napotkanego patternu kończacego

        if len(ending_places_right) == 0: # j.w. dot sekwencji "prawostronnej"
            end_place_right = len(right[n])
        else:
            end_place_right = ending_places_right[0]


        while i in range(end_place_left) and j in range(end_place_right): # funkcja bedzie działać dla zakresu i i j okreslonego powyżej

            if rev_left[i] == "(" and right[n][j] == ")":  # jeśli po obu stronach spinki nawiasy - identyfikacja jako pień

                if sequence_in_progress != 1:
                    struct[n].append({'type':'stem', 'start_l': i, 'start_r': j, 'length':1 })
                    sequence_in_progress = 1 # po wejściu do danej pętli np pnia, do listy dodawany jest nowy słownik zaw inf o strukturze
                elif sequence_in_progress == 1:
                    struct[n][-1]["length"] += 1 # przy kolejnym wejściu w tą samą pętlę dodawane jest jedynie +1 do wartości długości struktury
                i += 1
                j += 1 # przechodzenie symetryczne do kolejnych znaków


            elif rev_left[i] == "." and right[n][j] == ".": # jeśli po obu stronach spinki kropki - identyfikacja jako pętla

                if sequence_in_progress != 3:
                    struct[n].append({'type':'loop', 'start_l': i, 'start_r': j, 'length':0, 'assymetry_l': 0, 'assymetry_r':0 })
                    sequence_in_progress = 3 # jak poprzednio - dodanie nowego słownika do listy

                if sequence_in_progress == 3:

                    if rev_left[i+1] == "." and right[n][j+1] == ")":
                        struct[n][-1]["length"] += 1
                        struct[n][-1]["assymetry_l"] += 1 # jak poprzednio - dodanie +1 do dł struktury
                        i += 1 # w przypadku niesymetrycznej pętli dodawanie "i" i "j" niesymetryczne tj tylko tam gdzie "."
                    elif rev_left[i+1] == "(" and right[n][j+1] == ".":
                        struct[n][-1]["length"] += 1
                        struct[n][-1]["assymetry_r"] += 1
                        j += 1 # j. w.
                    else:
                        struct[n][-1]["length"] += 2 # działanie jak przy poprzedniej pętli
                        i += 1
                        j += 1 # jeśli pętla jest symetryczna - symetryczne dodawanie "i" i "j"


            elif rev_left[i] == "(" and right[n][j] == ".": # jeśli po jednej stronie spinki nawias a po 2giej kropka - identyfikacja jako bulge prawy

                if sequence_in_progress != 2: # jak poprzednio
                    struct[n].append({'type':'bulge', 'start_r': j, 'length':1 }) # bulge niesymetryczne, więc poczatek i koniec tylko z 1 strony
                    sequence_in_progress = 2
                elif sequence_in_progress == 2:
                    struct[n][-1]["length"] += 1
                j += 1 # bulge są niesymetryczne a więc niesymetryczne dodawanie "i " lub "j" tylko tam, gdzie "."

            elif rev_left[i] == "." and right[n][j] == ")": # jak poprzednio

                if sequence_in_progress != 2: # jak poprzednio
                    struct[n].append({'type':'bulge', 'start_l': i, 'length':1 })
                    sequence_in_progress = 2
                elif sequence_in_progress == 2:
                    struct[n][-1]["length"] += 1
                i += 1 # jak poprzednio

            else:
                break
    return struct

linear_structures = find_slb(left_structures_side, right_structures_side, structures)
print(linear_structures)


###################################################################################

# normalizacja numerów indeksów aby zgadzały się z indeksami wejściowej sekwencji

def index_normalization(left, right, linear_structures):
    for k in range(len(linear_structures)):
        for structure in linear_structures[k]:

            if "start_l" in structure.keys(): # warunek dodany, gdyż spinka ma tylko 1 początek i koniec a bulge mają długość tylko z 1 strony

                if structure['type'] == 'loop': # długość spinki jest liczona "podwójnie" dlatego konieczna jest pętla do poprawy tego
                    if structure['assymetry_l'] == 0 and structure['assymetry_r'] == 0:
                        structure['end_l'] = structure['start_l'] + (structure['length']/2)
                        # jeśli pętla symetryczna koniec z każdej strony to długośc/2
                    elif structure['assymetry_l'] != 0:
                        structure['end_l'] = structure['start_l'] + (structure['length']-((structure['length']- structure['assymetry_l'])/2))
                        # jeśli pętla symetryczna koniec to długośc/2
                    elif structure['assymetry_r'] != 0:
                        structure['end_l'] = structure['start_l'] + ((structure['length']- structure['assymetry_l'])/2)
                else:
                    structure['end_l'] = structure['start_l'] + structure['length']
                    # dopisanie do kazdego słownika (spinki/pętli/pnia/bulge) indeks zakończenia struktury po lewej stronie (start + length)

                structure['start_l'] = len(left[k]) - structure['start_l']
                structure['end_l'] = len(left[k]) - structure['end_l']
                # normalizacja indeksów sekwencji
                # sekwencje lewe były indeksowane na odwróconym stringu, dlatego od jego dł odejmuje się indeks

            if "start_r" in structure.keys(): # analogiczna pętla dla części prawych

                if structure['type'] == 'loop':
                    if structure['assymetry_l'] == 0 and structure['assymetry_r'] == 0:
                        structure['end_r'] = structure['start_r'] + (structure['length']/2)
                    elif structure['assymetry_r'] != 0:
                        structure['end_r'] = structure['start_r'] + (structure['length']-((structure['length']- structure['assymetry_l'])/2))
                    elif structure['assymetry_l'] != 0:
                        structure['end_r'] = structure['start_r'] + ((structure['length']- structure['assymetry_l'])/2)
                else:
                    structure['end_r'] = structure['start_r'] + structure['length']


                structure['start_r'] = structure['start_r'] + linear_structures[k][0]["end"]
                structure['end_r'] = structure['end_r'] + linear_structures[k][0]["end"]
                # normalizacja indeksów sekwencji
                # sekwencje prawe zaczynały się w raz z końcem danej spinki, do tej wartości trzeba dodac indeks

    return linear_structures


linear_structures = index_normalization(left_structures_side, right_structures_side, linear_structures)

print("Po znormalizowaniu indeksów;", linear_structures)

#############################################

# zamiana odpowiednich znaków w sekwencji na nazwy struktur : A0, A1 itp
# tylko dla struktur "1 rzędu" - znalezionych przy 1 iteracji
# przy kolejnych musza być dopisywane inne nazwy B0, B1 itd TODO

def add_linear_structure_names(linear_structures): # zapisuje wszystkie struktury wchodzące w skład str liniowej jako 1 obiekt
    linear_structures_names = []
    for p in range(len(linear_structures)):
        start_place = linear_structures[p][-1]["end_l"]
        end_place = linear_structures[p][-1]["end_r"]
        linear_structures_names.append({"name":"A"+str(p), "start":start_place, "end":end_place})
    # stworzenie nowej listy w której każda struktura liniowa otrzyma słownik z kluczami: imieniem, początkim i końcem
    return linear_structures_names


linear_structures_names = add_linear_structure_names(linear_structures)
print("Struktury liniowe jako nazwy i ich zakresy:", linear_structures_names)
print()

#############################################
# wizualizacja - zamiana odpowiednich znakow na nazwy struktur liniowych

def linear_structures_visualization(input_seq, linear_structures_names): # wizualizuje w jakich miejscach oryginalnej sekwencji znajdują się nazwy str liniowych
    dot_bracket_seq_2nd = ""

    for p in range(len(linear_structures_names)):
        if len(linear_structures_names) == 1:
            dot_bracket_seq_2nd += input_seq[:linear_structures_names[p]["start"]] + linear_structures_names[p]["name"] + input_seq[linear_structures_names[p]["end"]:]
            # jeśli tylko 1 struktura liniowa, nowa sekwencja do początku do startu struktury i od końca struktury do końca sekwencji
        else:
            if p == 0:
                dot_bracket_seq_2nd += input_seq[:linear_structures_names[p]["start"]] + linear_structures_names[p]["name"] + \
                                       input_seq[linear_structures_names[p]["end"]:linear_structures_names[p+1]["start"]]

            if p != 0 and p != (len(linear_structures_names)-1):
                dot_bracket_seq_2nd += linear_structures_names[p]["name"] + input_seq[linear_structures_names[p]["end"]:linear_structures_names[p+1]["start"]]

            if p == (len(linear_structures_names)-1):
                dot_bracket_seq_2nd += linear_structures_names[p]["name"] + input_seq[linear_structures_names[p]["end"]:]
                # początek sekwencji do startu 1szej (indeks 0) struktury liniowej, nazwa struktury, zakres od jej końca do kolejnej
                # (razy l struktur) i zakres do końca sekwencji
    return dot_bracket_seq_2nd


input_sequence_2nd = linear_structures_visualization(input_sequence, linear_structures_names)

print("Nowa wersja oryginalnej sekwencji:", input_sequence_2nd)

###########################################################

# odnajdywanie skrzyżowań - kropki między strukturami liniowymi

# znajdywanie znaków mdz strukturami
# TO NIE DZIAŁA DLA B SKOMPLIKOWANYCH TODO

if len(linear_structures_names) > 1: # ponieważ nie można wtedy znaleźć znaków pomiędzy
    junction_1_list = [] # tu będą zapisywane poszczególne skrzyżowania
    junction_in_progress = 0 # zmienna potrzebna do zaznaczenia wejścia do pętli

    for p in range(len(linear_structures_names)-1): # porównanie każdej struktury z kolejną
        print("to jest p", p)

        if linear_structures_names[p]['end'] != linear_structures_names[p+1]['start']: # gdy struktury są przedzielone jakimiś znakami
            print("znaki mdz strukturami")

            for sign in dot_bracket_seq[linear_structures_names[p]['end']: linear_structures_names[p+1]['start']]: # porównanie znaków mdz każdą str liniowa a kolejną
                print("To zakres pomiędzy:", linear_structures_names[p]['end'], ":", linear_structures_names[p+1]['start'] )
                print("to te znaki pomiędzy!", dot_bracket_seq[(linear_structures_names[p]['end']): (linear_structures_names[p+1]['start'])]) # TO JEST ZŁY ZNAK TODO
                print(junction_in_progress)

                if sign == ".":
                    print("tu jest kropka")

                    if junction_in_progress == 0: # przy wejściu do pętli dodawana jest nowa lista ze słownikiem ozn skrzyżowanie wraz z inf o długości i zawartości (nazwy str liniowych)

                        if junction_1_list == []: # jeśli jest to 1sza iteracja i jeszcze nie dodano elementów do listy
                            junction_1_list.append([{'name': 'J'+str(p), 'start': 0, 'end': 0, 'type':'junction', 'length': 1, 'content':[linear_structures_names[p], linear_structures_names[p+1]]}])

                        else:
                            if junction_1_list[-1][0]['type'] != 'junction': # jeśli dana struktura nie jest już zapisana w poprzednim skrzyżowaniu
                                junction_1_list.append([{'name': 'J'+str(p),'start': 0, 'end': 0, 'type':'junction', 'length': 1, 'content':[linear_structures_names[p], linear_structures_names[p+1]]}])
                                # utworzenie nowej listy zaw slownik ze skrzyżowaniem i inf o nim

                            elif junction_1_list[-1][0]['type'] == 'junction': # jeśli dana struktura jest już zapisana w poprzednim skrzyżowaniu
                                junction_1_list[-1][0]['length'] += 1 # jesli mdz str p(już zapisaną) a p+1 są same kropki - dodawane sa one do dł skrzyżowania
                                junction_1_list[-1][0]['content'].append(linear_structures_names[p+1]) # jesli mdz str p(już zapisaną) a p+1 są same kropki - dopisywana jest ona do zaw
                        junction_in_progress = 1

                    elif junction_in_progress == 1:
                        junction_1_list[-1][0]['length'] += 1

                else:
                    print("tu jest nawias")
                    if junction_1_list == []: # jeśli jest to 1sza iteracja i jeszcze nie dodano elementów do listy a 1 struktura nie tworzy skrzyżowania
                        junction_1_list.append([linear_structures_names[p]])

                    else:
                        if junction_1_list[-1][0]['type'] != 'junction': # jeśli dana struktura nie jest już zapisana w poprzednim skrzyżowaniu # A TU WYSKAKUJE BŁĄD TODO
                            junction_1_list.append([linear_structures_names[p]])

                        elif junction_1_list[-1][0]['type'] == 'junction': # jeśli dana struktura jest już zapisana w poprzednim skrzyżowaniu
                            junction_1_list.append(linear_structures_names[p+1])

                    break # jeśli między strukturami choć 1 znak inny niż "." to nie należą one do tego samego skrzyżowania
                print(junction_1_list)

        else:
            print("brak znakow mdz strukturami")
            if junction_1_list == []: # jeśli jest to 1sza iteracja i jeszcze nie dodano elementów do listy
                junction_1_list.append([{'name': 'J'+str(p),'type':'junction', 'start': 0, 'end': 0, 'length': 0, 'content':[linear_structures_names[p], linear_structures_names[p+1]]}])

            else:
                if junction_1_list[-1][0]['type'] != 'junction':
                    junction_1_list.append([{'name': 'J'+str(p),'type':'junction', 'start': 0, 'end': 0, 'length': 0, 'content':[linear_structures_names[p], linear_structures_names[p+1]]}])

                elif junction_1_list[-1][0]['type'] == 'junction':
                    junction_1_list[-1][0]['content'].append(linear_structures_names[p+1])

            print(junction_1_list) # jesli mdz strukturami nie ma znaków, dodawane są one do skrzyżowania, jednak bez zwiększenia jego długości

print("i ostatecznie:", junction_1_list)



# ##################################
#
# # dodawanie do skrzyżowania kropek znajdujących się po prawej i lewej stronie
# # BYĆ MOŻE MOGĄ BYĆ TEŻ KROPKI BEZPOŚREDNIO PRZYLEGAJĄCE DO SKRZYŻOWANIA
# # ALE NALEŻĄCE DO JAKIEJŚ INNEJ STRUKTURY - TRZEBA TO SPRAWDZIĆ TODO
# for p in range(len(junction_1_list)):
#     if junction_1_list[p][0]['type'] == 'junction':
#         junction_1_list[p][0]['start'] = junction_1_list[p][0]['content'][0]['start']
#         junction_1_list[p][0]['end'] = junction_1_list[p][0]['content'][-1]['end']
#         print("dodane start i end:", junction_1_list)
#
#         left_junction_side = dot_bracket_seq[:junction_1_list[p][0]['start']]
#         rev_left_junction_side = left_junction_side[::-1]
#         right_junction_side = dot_bracket_seq[junction_1_list[p][0]['end']:]
#         print("lewa:", left_junction_side)
#         print("lewa odwrócona:", rev_left_junction_side)
#         print("prawa:", right_junction_side)
#
#         for i in range (len(rev_left_junction_side)):
#             if rev_left_junction_side[i] == ".":
#                 print("to jest i:")
#                 junction_1_list[p][0]['length'] += 1
#                 junction_1_list[p][0]['start'] -= 1
#                 print("no i dodajemy z lewej")
#
#             else:
#                 print("i koniec z lewej")
#                 break
#
#         for j in range (len(right_junction_side)):
#             if right_junction_side[j] == ".":
#                junction_1_list[p][0]['length'] += 1
#                junction_1_list[p][0]['end'] += 1
#                print("no i dodajemy z prawej")
#             else:
#                 print("i koniec z lewej")
#                 break
#
# print("i po dodaniu kropek obok:", junction_1_list)


#########################################################
# zapisanie wszystkich dotychczasowych rzeczy jako funkcji # TODO
# zrobienie dużej pętli z funkcjami pracującej aż do końca # TODO
# odpowiednie zakończenie pętli # TODO

# do sekwencji wyzejrzędowych nie dopisywac A1, A2 itp, tylko całe struktury zagnieżdżone w sobie  (???) TODO

# odszukiwanie indeksów pseudowęzłów i przypisywanie ich do konkretnych struktur TODO

# deal with igraph TODO



####################################################################3
# print()
# print("proba dodawania")
# lista_list = [[{1:'a', 2:'b'}, {3:'c', 4:'d'}], [{'z':24, 'y':23}, {'x':22, 'w':21}], ["hey", "you"]]
# lista_list[0].append({'f':17})
# print("Lista", lista_list)
# lista_list[1][0]["z"] += 1
# print("Lista2", lista_list)
# print(lista_list[2][0][1])
#
# def probna_funkcja(lista):
#     if lista[1][0]["y"] == 23:
#         print("superowo!")
#     else:
#         print("nie do końca")
#
# probna_funkcja(lista_list)
#
# # print(len(left_hairpin_side[-1]))
#
# # for i in range(0,100):
# #     if i in range(0,5) and i in range(0,3):
# #         print(i)
#
# k = [1,4,5,6,7,8]
# p = []
# print(len(p))

# pustalista = []
# for f in range(5):
#     print(f)
#     pustalista.append("A"+str(f))
# print("Test", pustalista)

# listka = [0]
# print(len(listka))
# for p in range(len(listka)):
#     print(p)
# print(len(listka)-1)

