# -*- coding: utf-8 -*-
# zadanie:
# zamiana struktury 2nd RNA zapisane w formacie dot-bracket na graf przedstawiający rozmieszczenie poszczególnych struktur

import re
import igraph

def find_seq(input_file): # zwraca sekwencję w formacie dot-bracket z wczytywanego pliku
    sequence = ""
    bases = ["A", "C", "U", "G", "N"]
    signs = [".", "(", ")", "[", "]", ">", "<", "{", "}"]
    # listy zasad azotowych (po nich następuje właściwa sekwencja dot-bracket), oraz znaków w poszukiewanym formacie
    for i in range(len(input_file)):
        if input_file[i] in bases and input_file[i+1] in signs and input_file[i+2] in signs and input_file[i+3] in signs:
            sequence = (input_file[i+1:])
    return sequence

def collecting_input(): #otwiera wczytywany plik lub sekwencję i przekształca w sekwencję wejściową dla programu (w formacie dot-bracket)
    print("If you want to load a file from RNA Soft STRAND database: print 'file'")
    print("If you want to paste sequence in dot-bracket (dot-parenthesis) format: print 'sequence'")
    file_or_seq = raw_input("> ")

    if file_or_seq == 'file':
        print("Please paste path to the file")
        file_path = raw_input("> ")
        in_file = open(file_path)
        in_file2 = in_file.read().replace("\n", "")
        input_file = str(in_file2) # importowanie pliku wejściowego, usunięcie newline i zamiana na string
        dot_bracket_seq = find_seq(input_file)

    elif file_or_seq == 'sequence':
        print("Please paste the sequence in dot-bracket format")
        input_seq = raw_input("> ")
        dot_bracket_seq = input_seq.replace("\n", "") # NIE DZIAŁA KIEDY JEST SEKWENCA W WIELU LINIJKACH TODO

    else:
        print("Please choose between 'file' and 'sequence'")
        collecting_input()

    return dot_bracket_seq

dot_bracket_seq = collecting_input()

print("Sekwencja w formacie dot-bracket:", dot_bracket_seq)


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


hairpin_structures = find_hairpins(input_sequence)
print("Znalezione struktury spinek do włosów:", hairpin_structures)
print()

################################################################################


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

left_structures_side, right_structures_side = left_and_right_structure_side(input_sequence, hairpin_structures)

print("Lista sekwencji po lewej stronie struktur:", left_structures_side)
print("Lista sekwencji po prawej stronie struktur:", right_structures_side)
print()

#################################################################################################


def find_slb(left, right, structures_list):  # odnajduje pnie/pętle/bulge wew(stems/bulges/loops) połączone liniowo ze spinkami lub str wyższego rzędu
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
                    structures_list[n].append({'type':'stem', 'start_l': i, 'start_r': j, 'length':1, 'start':0, 'end':0 })
                    sequence_in_progress = 1 # po wejściu do danej pętli np pnia, do listy dodawany jest nowy słownik zaw inf o strukturze
                elif sequence_in_progress == 1:
                    structures_list[n][-1]["length"] += 1 # przy kolejnym wejściu w tą samą pętlę dodawane jest jedynie +1 do wartości długości struktury
                i += 1
                j += 1 # przechodzenie symetryczne do kolejnych znaków


            elif rev_left[i] == "." and right[n][j] == ".": # jeśli po obu stronach spinki kropki - identyfikacja jako pętla

                if sequence_in_progress != 3:
                    structures_list[n].append({'type':'loop', 'start_l': i, 'start_r': j, 'length':0, 'assymetry_l': 0, 'assymetry_r':0, 'start':0, 'end':0 })
                    sequence_in_progress = 3 # jak poprzednio - dodanie nowego słownika do listy

                if sequence_in_progress == 3:

                    if rev_left[i+1] == "." and right[n][j+1] == ")":
                        structures_list[n][-1]["length"] += 1
                        structures_list[n][-1]["assymetry_l"] += 1 # jak poprzednio - dodanie +1 do dł struktury
                        i += 1 # w przypadku niesymetrycznej pętli dodawanie "i" i "j" niesymetryczne tj tylko tam gdzie "."
                    elif rev_left[i+1] == "(" and right[n][j+1] == ".":
                        structures_list[n][-1]["length"] += 1
                        structures_list[n][-1]["assymetry_r"] += 1
                        j += 1 # j. w.
                    else:
                        structures_list[n][-1]["length"] += 2 # działanie jak przy poprzedniej pętli
                        i += 1
                        j += 1 # jeśli pętla jest symetryczna - symetryczne dodawanie "i" i "j"


            elif rev_left[i] == "(" and right[n][j] == ".": # jeśli po jednej stronie spinki nawias a po 2giej kropka - identyfikacja jako bulge prawy

                if sequence_in_progress != 2: # jak poprzednio
                    structures_list[n].append({'type':'bulge', 'start_r': j, 'length':1, 'start':0, 'end':0 }) # bulge niesymetryczne, więc poczatek i koniec tylko z 1 strony
                    sequence_in_progress = 2
                elif sequence_in_progress == 2:
                    structures_list[n][-1]["length"] += 1
                j += 1 # bulge są niesymetryczne a więc niesymetryczne dodawanie "i " lub "j" tylko tam, gdzie "."

            elif rev_left[i] == "." and right[n][j] == ")": # jak poprzednio

                if sequence_in_progress != 2: # jak poprzednio
                    structures_list[n].append({'type':'bulge', 'start_l': i, 'length':1, 'start':0, 'end':0 })
                    sequence_in_progress = 2
                elif sequence_in_progress == 2:
                    structures_list[n][-1]["length"] += 1
                i += 1 # jak poprzednio

            else:
                break
    return structures_list

linear_structures = find_slb(left_structures_side, right_structures_side, hairpin_structures)
print(linear_structures)


###################################################################################

def index_normalization(left, right, linear_structures): # normalizuje numery indeksów aby zgadzały się z indeksami wejściowej sekwencji
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
                    structure['start'] = structure['end_l']
                elif structure['type'] == 'stem':
                    structure['end_l'] = structure['start_l'] + structure['length']
                    structure['start'] = structure['end_l']
                else:
                    structure['end_l'] = structure['start_l'] + structure['length']
                    structure['start'] = structure['end_l']
                    structure['end'] = structure['start_l']
                    structure['end'] = len(left[k]) - structure['end']
                    # dopisanie do kazdego słownika (spinki/pętli/pnia/bulge) indeks zakończenia struktury po lewej stronie (start + length)

                structure['start_l'] = len(left[k]) - structure['start_l']
                structure['end_l'] = len(left[k]) - structure['end_l']
                structure['start'] = len(left[k]) - structure['start']
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
                    structure['end'] = structure['end_r']
                elif structure['type'] == 'stem':
                    structure['end_r'] = structure['start_r'] + structure['length']
                    structure['end'] = structure['end_r']
                else:
                    structure['end_r'] = structure['start_r'] + structure['length']
                    structure['start'] = structure['start_r']
                    structure['end'] = structure['end_r']
                    structure['start'] = structure['start'] + linear_structures[k][0]["end"]

                structure['start_r'] = structure['start_r'] + linear_structures[k][0]["end"]
                structure['end_r'] = structure['end_r'] + linear_structures[k][0]["end"]
                structure['end'] = structure['end'] + linear_structures[k][0]["end"]
                # normalizacja indeksów sekwencji
                # sekwencje prawe zaczynały się w raz z końcem danej spinki, do tej wartości trzeba dodac indeks

    return linear_structures

linear_structures = index_normalization(left_structures_side, right_structures_side, linear_structures)

print("Struktury liniowe (po znormalizowaniu indeksów):", linear_structures)

#############################################

def add_linear_structures_names(linear_structures): # zapisuje wszystkie struktury wchodzące w skład str liniowej jako 1 obiekt
    linear_structures_names = []
    for p in range(len(linear_structures)):
        start_place = linear_structures[p][-1]["end_l"]
        end_place = linear_structures[p][-1]["end_r"]
        linear_structures_names.append({"name":"A"+str(p), "start":start_place, "end":end_place})
    # stworzenie nowej listy w której każda struktura liniowa otrzyma słownik z kluczami: imieniem, początkim i końcem
    return linear_structures_names

linear_structures_names = add_linear_structures_names(linear_structures)
print("Struktury liniowe jako nazwy i ich zakresy:", linear_structures_names)
print()

#############################################


def linear_structures_visualization(input_seq, linear_structures_names): # wizualizuje w jakich miejscach oryginalnej sekwencji znajdują się nazwy str liniowych
    linear_str_vis = ""

    for p in range(len(linear_structures_names)):
        if len(linear_structures_names) == 1:
            linear_str_vis += input_seq[:linear_structures_names[p]["start"]] + linear_structures_names[p]["name"] + input_seq[linear_structures_names[p]["end"]:]
            # jeśli tylko 1 struktura liniowa, nowa sekwencja do początku do startu struktury i od końca struktury do końca sekwencji
        else:
            if p == 0:
                linear_str_vis += input_seq[:linear_structures_names[p]["start"]] + linear_structures_names[p]["name"] + \
                                       input_seq[linear_structures_names[p]["end"]:linear_structures_names[p+1]["start"]]

            if p != 0 and p != (len(linear_structures_names)-1):
                linear_str_vis += linear_structures_names[p]["name"] + input_seq[linear_structures_names[p]["end"]:linear_structures_names[p+1]["start"]]

            if p == (len(linear_structures_names)-1):
                linear_str_vis += linear_structures_names[p]["name"] + input_seq[linear_structures_names[p]["end"]:]
                # początek sekwencji do startu 1szej (indeks 0) struktury liniowej, nazwa struktury, zakres od jej końca do kolejnej
                # (razy l struktur) i zakres do końca sekwencji
    return linear_str_vis

linear_structures_visual = linear_structures_visualization(input_sequence, linear_structures_names)

print("Wizualizacja struktur liniowych:", linear_structures_visual)

###########################################################

def junction_find(input_sequence, linear_structures): # odnajduje struktury liniowe połączone bezpośrednio lub "." - będące w 1 skrzyżowaniu
    junction_list = [] # tu będą zapisywane poszczególne skrzyżowania

    if len(linear_structures) > 1: # ponieważ nie można wtedy porównać kolejnych struktur linowych a pojedyncze skrzyżowanie nie tworzy wyrzejrzędowych skrzyżowań
        for p in range(len(linear_structures)-1): # porównanie każdej struktury z kolejną aż do przedostatniej

            if linear_structures[p][-1]['end'] != linear_structures[p+1][-1]['start']: # gdy struktury są przedzielone jakimiś znakami TUTAJ -1 i -1

                if all(sign == "." for sign in input_sequence[linear_structures[p][-1]['end']: linear_structures[p+1][-1]['start']] ) == True: #TUTAJ -1 i -1
                    # kiedy wszystkie znaki pomiędzy strukturami są kropkami

                    if junction_list == []:
                        junction_list.append([{'name': 'J'+str(p), 'start': 0, 'end': 0, 'type':'junction',
                                                 'length': (linear_structures[p+1][-1]['start'] - linear_structures[p][-1]['end']) , # TUTAJ -1 i -1
                                                 'content':[linear_structures[p], linear_structures[p+1]]}])
                        # jeśli jest to 1sza iteracja i jeszcze nie dodano elementów do listy a pierwsze 2 struktury tworzą skrzyżowanie
                        # dodawana jest lista ze słownikiem z nowym skrzyżowaniem zaw struktury p i p+1
                    else:
                        junction_list[-1][0]['content'].append(linear_structures[p+1])
                        junction_list[-1][0]['length'] += (linear_structures[p+1][-1]['start'] - linear_structures[p][-1]['end']) # TUTAJ -1 i -1
                        # jeśli kolejna znaleziona str należy do skrzyzowania dopisywana jest do poprzedniego wraz z ilością kropek pomiędzy

                else: # kiedy między strukturami jest min i ")" lub ")" - struktury p i p+1 nie sa w 1 skrzyżowaniu
                    if junction_list == []:
                        junction_list.append([{'name': 'J'+str(p), 'start': 0, 'end': 0, 'type':'junction',
                                                 'length': 0, 'content':[linear_structures[p]]}])
                        junction_list.append([{'name': 'J'+str(p+1), 'start': 0, 'end': 0, 'type':'junction',
                                                 'length': 0, 'content':[linear_structures[p+1]]}])
                        # jeśli jest to 1sza iteracja i jeszcze nie dodano elementów do listy a pierwsze 2 struktury nie są w 1 skrzyżowaniu
                        # dodawane są 2 nowe listy ze słownikami ze skrzyżowaniami zawierającymi 1) - strukturę p i 2) - strukturę p+1
                    else:
                        junction_list.append([{'name': 'J'+str(p+1), 'start': 0, 'end': 0, 'type':'junction',
                                                 'length': 0, 'content':[linear_structures[p+1]]}])
                        # jeśli kolejna str nie należy do poprzedniego skrzyżowania zapisywana jest jako nowe zaw str p+1

            else: # jeśli mdz strukturami nie ma żadnego znaku - należą one do 1 skrzyżowania
                if junction_list == []:
                    junction_list.append([{'name': 'J'+str(p), 'start': 0, 'end': 0, 'type':'junction',
                                             'length': 0, 'content':[linear_structures[p], linear_structures[p+1]]}])
                else:
                    junction_list[-1][0]['content'].append(linear_structures[p+1])
                   # jak w pętli dla znalezienia ".", poza dodawaniem długości

    return junction_list


junctions = junction_find(input_sequence, linear_structures)

print("To są początki skrzyżowań)", junctions)

##################################


def junction_elongation(input_sequence, junction_list): # dodaje do skrzyżowania kropki znajdujące się po prawej i lewej stronie
    for p in range(len(junction_list)):
        junction_list[p][0]['start'] = junction_list[p][0]['content'][0][-1]['start'] # początek skrzyżowania to początek jego 1szej struktury liniowej TUTAJ -1
        junction_list[p][0]['end'] = junction_list[p][0]['content'][-1][-1]['end'] # koniec skrzyżowania to koniec jego ostatniej struktury liniowej TUTAJ -1

        left_junction_side = input_sequence[:junction_list[p][0]['start']] # sekwencja po lewej stronie od skrzyżowania
        rev_left_junction_side = left_junction_side[::-1] # odwrócona "lewa" sekwencja - aby łatwiejsza była iteracja
        right_junction_side = input_sequence[junction_list[p][0]['end']:] # sekwencja po prawej stronie od skrzyżowania

        for i in range (len(rev_left_junction_side)):
            if rev_left_junction_side[i] == ".": # jeśli znak po lewej stronie od skrzyżowania to "."
                junction_list[p][0]['length'] += 1 # dodawane jest +1 do długości skrzyżowania
                junction_list[p][0]['start'] -= 1 # początek skrzyżowania jest przesuwany o 1 w lewą stronę

            else: # jeśli funkcja napotka na nawias - koniec skrzyżowania
                break

        for j in range (len(right_junction_side)): # analogiczna pętla dla prawej strony od skrzyżowania
            if right_junction_side[j] == ".":
               junction_list[p][0]['length'] += 1
               junction_list[p][0]['end'] += 1
            else:
                break
    return junction_list


junctions = junction_elongation(input_sequence, junctions)
print("Odnalezione skrzyżowania:", junctions)
print("Jest ich", len(junctions))
print()

#########################################################################
# wizualizacja skrzyżowań

def junction_visualization(input_seq, junction_list): # wizualizuje w jakich miejscach oryginalnej sekwencji znajdują się nazwy str liniowych
    junction_vis = ""

    for p in range(len(junction_list)):
        if len(junction_list) == 1:
            junction_vis += input_seq[:junction_list[p][0]["start"]] + junction_list[p][0]["name"] + input_seq[junction_list[p][0]["end"]:]
            # jeśli tylko 1 skrzyżowanie, nowa sekwencja do początku do startu skrzyżowania i od końca skrzyżowania do końca sekwencji
        else:
            if p == 0:
                junction_vis += input_seq[:junction_list[p][0]["start"]] + junction_list[p][0]["name"] + \
                                       input_seq[junction_list[p][0]["end"]:junction_list[p+1][0]["start"]]

            if p != 0 and p != (len(junction_list)-1):
                junction_vis += junction_list[p][0]["name"] + input_seq[junction_list[p][0]["end"]:junction_list[p+1][0]["start"]]

            if p == (len(junction_list)-1):
                junction_vis += junction_list[p][0]["name"] + input_seq[junction_list[p][0]["end"]:]
                # początek sekwencji do startu 1szej (indeks 0) struktury liniowej, nazwa struktury, zakres od jej końca do kolejnej
                # (razy liczba struktur) i zakres do końca sekwencji
    return junction_vis


junctions_visual = junction_visualization(input_sequence, junctions)
print("Wizualizacja skrzyżowań:", junctions_visual)

########################################################
########################################################
# funkcje potrzebne w szukaniu struktur i skrzyżowań 2go i dalszorzedowych
# odnajdywanie skrzyzowań na końcach struktur liniowych które następnie będą rozbudowane do struktury liniowej
def find_hairpin_junctions(input_seq, left, right, junctions): # funkcja odnajdująca skrzyżowania będące odpowiednikami spinek do włosów tj na końcu str liniowych
    hairpin_junctions_list = []
    hairpin_junctions_left = []
    hairpin_junctions_right = []
    for n in range(len(junctions)):
        if junctions[n][0]['start'] != 0 or junctions[n][0]['end'] != len(input_seq): # skrzyżowanie na początku i końcu sekwencji nie może być odpowiednikiem spinki
            if input_seq[junctions[n][0]['start']-1] == "(" and input_seq[junctions[n][0]['end']+1] == ")": # skrzyżowanie-spinka jest otoczone "(" po lewej i ")" po prawej str
                hairpin_junctions_list.append(junctions[n])
                hairpin_junctions_left.append(left[n])
                hairpin_junctions_right.append(right[n])
    return hairpin_junctions_list, hairpin_junctions_left, hairpin_junctions_right # tworzenie list ze skrzyżowaniami-sponkami i ic lewymi i prawymi stronami

##################
# lista wszystkich skrzyżowań - tych rozbudowywanych i tych nie
def new_structures_list(junctions, linear_structures): #tworzenie listy z pozostałymi skrzyżowaniami i nowymi strukturami liniowymi
    for structure in junctions:
        for n in range(len(linear_structures)):
            if structure[0]['name'] == linear_structures[n][0]['name']:
                structure == linear_structures[n]
    return junctions

#######################################################
#######################################################
# funkcja szukająca 2go i kolejnorzędowych struktur liniowych i skrzyżowań aż do końca sekwencji
iteration_number = 1
def the_ultimate_function(input_seq, junctions, iteration_number):
    # iteration_number += 1

    if junctions[0][0]['start'] == 0 and junctions[0][0]['end'] == len(input_seq):
        print("koniec struktury")
        return junctions, iteration_number
    else:
        iteration_number += 1
        print()
        print("szukanie struktur i skrzyżowań", iteration_number, "rzędowych")
        left_junctions_side, right_junctions_side = left_and_right_structure_side(input_sequence, junctions)

        hairpin_junction_structures, hairpin_junction_left, hairpin_junction_right = find_hairpin_junctions(input_sequence, left_junctions_side, right_junctions_side, junctions)
        print("Skrzyżowania - spinki:", hairpin_junction_structures)

        linear_structures = find_slb(hairpin_junction_left, hairpin_junction_right, hairpin_junction_structures)
        linear_structures = index_normalization(hairpin_junction_left, hairpin_junction_right, linear_structures)
        linear_structures = new_structures_list(junctions, linear_structures)
        print(iteration_number, "rzędowe struktury liniowe:", linear_structures)

        if len(linear_structures) > 1: # tylko jesli jest >1 skrzyżowanie będzie ono dalej łączone w skrzyżowania wyżejrzędowe

            junctions = junction_find(input_sequence, linear_structures)
            junctions = junction_elongation(input_sequence, junctions)
            print(iteration_number, "rzędowe skrzyżowania:", junctions)
            print("Jest ich", len(junctions))

            junctions_visual = junction_visualization(input_sequence, junctions)
            print("Wizualizacja skrzyżowań", iteration_number, "rzędowych:", junctions_visual)

            return the_ultimate_function(input_seq, junctions, iteration_number)

        else:
            return linear_structures, iteration_number



# TODO nie działa dla example_mega : iteruje się w nieskonczoność...
# przy szukaniu struktur 3ciorzędowych nie tworzy się lista ze skrzyżowaniami-spinkami
# skrzyżowanie J0 ma strukturę liniową tylko po prawej stronie - nowa, dziwna struktura która nie jest ujęta w algorytmie - "spinka" na skrzyżowaniu 2ch spinek


final_list, iteration_number = the_ultimate_function(input_sequence, junctions, iteration_number)
print("Ostateczna lista skrzyżowań", final_list)
print("iteration number:", iteration_number)

###################################################################################
###################################################################################
# tworzenie grafu wizualizującego struktury liniowe
print()
print("WIZUALIZACJA")

def find_structures(junctions): # funkcja odnajdujaca struktury i zapisująca ich nazwy oraz długości do list
    v_names = [] #lista z nazwami struktur
    v_length = [] # lista z długościami struktur

    for structure in junctions:
        for n in range(len(structure)-1, -1, -1): #iteracja od tyłu, aby spinki znajdowały się "na końcu"

            if structure[n]['type'] == 'junction':
                v_names.append(structure[n]['type'])
                v_length.append(structure[n]['length'])

                junctions = structure[n]['content'] # jesli struktura jest skrzyżowaniem powtórne wywołanie funkcji na nowych argumentach

                temporary_v_names, temporary_v_length = find_structures(junctions)
                v_names.extend(temporary_v_names) # dodanie nazwy i długości skrzyżowania do odpowiednich list
                v_length.extend(temporary_v_length)

            else:
                v_names.append(structure[n]['type']) # dodanie nazwy i długości struktury do odpowiednich list
                v_length.append(structure[n]['length'])

    return v_names, v_length

v_names_raw, v_lengths_raw = find_structures(final_list)

#################################################
# usuwanie z listy pseudoskrzżowań - skrzyżowań bezposrednio połączonych ze sobą, które są artefaktem działania algorytmu
v_names = [] # nowa
v_lengths = []
for i in range(len(v_names_raw)):
    if v_names_raw[i] == "junction" and v_names_raw[i+1] == "junction": # jeśli 2 skrzyżowania są bezpośrednio obok siebie, 1 z nich nie zostaje dopisane do kolejnej listy
        None
    else:
        v_names.append(v_names_raw[i])
        v_lengths.append(v_lengths_raw[i])

print()
print("nazwy struktur", v_names)
print("długosci struktur", v_lengths)
print("ilość struktur", (len(v_names)))

#########################
# funkcja wyszukująca połączenia miedzy strukturami i dodająca krawędzie do listy
def find_edges(v_names):
    edges = []
    for n in range(len(v_names)):
        if v_names[n] == "stem": # pień jest bezpośrednio połączony liniowo ze strukturami po obu jego stronach

            if n != 0:
                if v_names[n-1] != "hairpin":
                    edges.append((n-1, n))

                else:
                    for k in range(0, n):
                        if v_names[k] == 'junction':
                            junction_for_stem = k

                    edges.append((junction_for_stem, n))

            edges.append((n, n+1))

        if v_names[n] == "junction": # skrzyżowanie jest połączone liniowo z kolejnym najbliższym z lewej skrzyżowaniem
            if n != 0:
                for k in range(0, n):
                    if v_names[k] == 'junction':
                        junction_for_stem = k

                edges.append((junction_for_stem, n))


    return edges


edges_list = find_edges(v_names)
print("lista krawędzi:", edges_list)

###############################################
# twrzenie grafu za pomocą pakietu igraph i wizualizacja za pomoca pakietu cairo

rna_graph = igraph.Graph(edges_list) # tworzenie grafu w oparciu o stworzoną listę krawędzi

rna_graph.vs["name"] = v_names # właściwości wierzchołków - nazwy struktur
rna_graph.vs["length"] = v_lengths # i ich długości

rna_graph_layout = rna_graph.layout("kk") # wybór wizualizacji (graf nieskierowany Kamada-Kawai)
igraph.plot(rna_graph, layout = rna_graph_layout) # rysowanie grafu


color_dict = {"stem":"brown", "hairpin":"red", "bulge":"green", "loop":"blue", "junction":"yellow"} # dodanie kodu kolorów oznaczających różne struktury

visual_style = {} # dodanie dodatkowych parametrów wizualnych grafu
visual_style["vertex_size"] = rna_graph.vs["length"] # wielkość wierzchołka zalezy od długości struktury
visual_style["vertex_color"] = [color_dict[name] for name in rna_graph.vs["name"]] #kolor wierzchołka zalezy od rodzaju struktury
visual_style["vertex_label"] = rna_graph.vs["name"] #metka wierzchołka jest rodzajem jego struktury
visual_style["layout"] = rna_graph_layout # wyspecyfikowany wcześniej rodzaj wizualizacji
igraph.plot(rna_graph, **visual_style) #graf zawierający wszystkie dodatkowe atrybuty



# odszukiwanie indeksów pseudowęzłów i przypisywanie ich do konkretnych struktur TODO


