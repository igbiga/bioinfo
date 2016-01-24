# zadanie:
# zamiana struktury 2nd RNA zapisane w formacie dot-bracket na graf przedstawiający rozmieszczenie poszczególnych struktur

import re
in_file = open("C:/Users/Igusia/Documents/example.txt")
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
print("Sekwencja w formacie dot-bracket: \n", dot_bracket_seq)
# wypisanie sekwencji w formacie d-b

#####################################
# 1 czesc algorytmu: wyszukiwanie struktur liniowych złozonch ze spinek, pni/bulge/pętli wewnętrznych


# wyszukiwanie spinek
# def find_hairpins(seq):
hairpin_pat = "\(\.*\)" # wyrazenie regularne znajdujace spinki (dowolna ilość "." zakończona "(" z lewej i ")" z prawej strony)
hairpin_pattern = re.compile(hairpin_pat) # kompilacja wyrażenia regularnego

hairpins = []
for m in hairpin_pattern.finditer(dot_bracket_seq): # pętla dla każdego patternu wyszukanego w sekwencji
    hairpins.append({'type': 'hairpin', 'start': m.start() + 1, 'end': m.end() - 1, 'length': m.end() - m.start()})
print("Znalezione struktury spinek:\n", hairpins)
# stworzenie listy słowników z których każdy będzie zawierał inf o kolejnych spinkach: rozpoczęcie, zakończenie i długość

########

structures = []
for n in range(len(hairpins)):
    structures.append([])
    #structures[n].append("L"+ n)
    structures[n].append(hairpins[n])
print(structures) # lista z poszczególnymi strukturami liniowymi.
# każda str jest zapisana w osobnej liście skł się ze słowników z informacjami
# o jej składnikach tj spinkach/pniach/pętlach/bulge, miejscach ich rozpoczecia i zakończenia w str i długości

#######

# wyszukiwanie pnia/bulge/pętli wew

left_hairpin_side = []
right_hairpin_side = []
for n in range(0, len(hairpins)):
    print("Spinka ",n, ":", hairpins[n])
    left_side = dot_bracket_seq[0: hairpins[n]['start']]
    left_hairpin_side.append(left_side) # lista sekwencji po lewej stronie każdej ze znalezionch spinek
    right_side = dot_bracket_seq[hairpins[n]['end']: -1]
    right_hairpin_side.append(right_side) # analogiczna lista z sekwencjami na prawo od spinek
# dla każdej znalezionej spinki wypisanie fragmentu sekwencji po jej lewej oraz prawej stronie
# bez uwzględnienia znaków wchodzących w skłąd spinki
# stworzenie list sekwencji "prawostronnych" oraz "lewostronnnych"

print("Lista sekwencji po lewej stronie spinek:", left_hairpin_side)
print("Lista sekwencji po prawej stronie spinek:", right_hairpin_side)

#################################################################################################
# funkcja odnajdująca pnie/pętle/bulge

def find_slb(left, right, struct):  # funkcja odnajdująca pnie/pętle/bulge wew(stems/bulges/loops) połączone liniowo ze spinkami
    #print("Lewa:", left, "i prawa:", right)

    ending_left_pat = "\.*\)" # wyrażenie regularne odnajdujące koniec lewostronny danej struktury liniowej (tj kropki i nawias w odwrotą stronę)
    ending_left_pattern = re.compile(ending_left_pat) # kompilacja wyrażenia regularnego
    ending_right_pat = "\.*\(" # jw. - koniec prawostronny
    ending_right_pattern = re.compile(ending_right_pat) # kompilacja wyrażenia regularnego


    for n in range(len(left)):  # zakres odp ilości "lewych części spinek" a więc i znalezionych spinek
        rev_left = left[n][::-1] # odwrócenie sekwencji "lewych cześci spinek" aby łatwiejsze było porównywanie sekwencji "prawej" i "lewej"
        #print("Odwrotność:", rev_left)
        #print("Dla ", n, "spinki część lewostronna ma długość: ", len(left[n]), " a prawostronna ma długość: ", len(right[n]))
        sequence_in_progress = 0  # zmienna oznaczająca strukturę (pień/pętlę/bulge) w której pętli obecnie znajduje się funkcja
        i = 0 # zmienna oznaczająca znak w sekwencji "lewej"
        j = 0 # zmienna oznaczająca znak w sekwencji "prawej"

        ending_places_left = [] # lista oznaczająca znalezione na podstawie patternu miejsca kończące funkcji po lewej stronie
        ending_places_right = [] # j.w - po prawej stronie

        for end_l in ending_left_pattern.finditer(rev_left): # pętla dla każdego patternu "kończącego" wyszukanego w sekwencji "lewostronnej"
            ending_places_left.append(end_l.start()) # dodanie indeksu znaku rozpoczynającego pattern do listy
            #print("Znalezione miejsca kończące w lewej części od spinki:\n", ending_places_left)
        for end_r in ending_right_pattern.finditer(right[n]): # j.w - dotyczy strony prawej
            ending_places_right.append(end_r.start())
            #print("Znalezione miejsca kończące w prawej części od spinki:\n", ending_places_right)


        if len(ending_places_left) == 0: # jeśli nie znaleziono żadnego patternu "konczącego" funkcja będzie działać aż do końca sekwencji
            end_place_left = len(rev_left)
        else:
            end_place_left = ending_places_left[0] # jeśli zostanie znaleziony pattern "kończący" funkcja będzie działać aż do 1szego napotkanego

        if len(ending_places_right) == 0: # j.w. dot sekwencji "prawostronnej"
            end_place_right = len(right[n])
        else:
            end_place_right = ending_places_right[0]

        #print("KONIEC LEWY", end_place_left, "I PRAWY", end_place_right)


        while i in range(end_place_left) and j in range(end_place_right): # funkcja bedzie działać dla zakresu i i j okreslonego powyżej
            # print("To jest i:", i, "a to j:", j)
            structure_length = len(struct[n])  # zmienna oznaczająca aktualną długość listy odpowiadającej strukturze liniowej do której dopisywane są kolejne el: pnie/pętle itp

            if rev_left[i] == "(" and right[n][j] == ")":  # jeśli po obu stronach spinki nawiasy - identyfikacja jako pień
                # print(rev_left[i], "oraz",  right[n][j], "to pien")

                if sequence_in_progress != 1:
                    struct[n].append({'type':'stem', 'start_l': i, 'start_r': j, 'length':1 })
                    sequence_in_progress = 1 # po wejściu do danej pętli np pnia, do listy dodawany jest nowy słownik zaw inf o strukturze
                elif sequence_in_progress == 1:
                    struct[n][-1]["length"] += 1 # przy kolejnym wejściu w tą samą pętlę dodawane jest jedynie +1 do wartości długości
                i += 1
                j += 1 # przechodzenie symetryczne do kolejnych znaków
                # print("calutka lista:", struct)
                # print ("Nowe i:", i , "oraz j:", j)
                # print("dlugosc struktury", structure_length)

            elif rev_left[i] == "." and right[n][j] == ".": # jeśli po obu stronach spinki kropki - identyfikacja jako pętla
                # print(rev_left[i], "oraz",  right[n][j], "to petla")

                if sequence_in_progress != 3:
                    struct[n].append({'type':'loop', 'start_l': i, 'start_r': j, 'length':2, 'assymetry': 'no' })
                    sequence_in_progress = 3
                    i += 1
                    j += 1
                elif sequence_in_progress == 3:

                    if rev_left[i+1] == "(" and right[n][j+1] == ".":
                        # print(rev_left[i+1], "oraz",  right[n][j+1], "to petla ale niesymetryczna z lewej")
                        struct[n][-1]["length"] += 1
                        struct[n][-1]["assymetry"] = 'right'
                        j += 1 # w przypadku niesyetrycznej pętli dodawanie "i" i "j" niesymetryczne tj tylko tam gdzie "."
                    elif rev_left[i+1] == "." and right[n][j+1] == ")":
                        # print(rev_left[i+1], "oraz",  right[n][j+1], "to petla ale niesymetryczna z prawej")
                        struct[n][-1]["length"] += 1
                        struct[n][-1]["assymetry"] = 'left'
                        i += 1 # j.w
                    else:
                        struct[n][-1]["length"] += 2 # działanie jak przy poprzedniej pętli
                        i += 1
                        j += 1 # jeśli pętla jest symetryczna - symetryczne dodawanie "i" i "j"


            elif rev_left[i] == "(" and right[n][j] == ".": # jeśli po jednej stronie spinki nawias a po 2giej kropka - identyfikacja jako bulge
                # print(rev_left[i], "oraz",  right[n][j], "to bulka prawa")

                if sequence_in_progress != 2: # jak poprzednio
                    struct[n].append({'type':'bulge', 'start_r': j, 'length':1 })
                    sequence_in_progress = 2
                elif sequence_in_progress == 2:
                    struct[n][-1]["length"] += 1
                j += 1 # bulge są niesymetryczne a więc niesymetryczne dodawanie "i " lub "j" tylko tam, gdzie "."
                # print("calutka lista:", struct)
                # print("dlugosc struktury", structure_length)

            elif rev_left[i] == "." and right[n][j] == ")": # jak poprzednio
                # print(rev_left[i], "oraz",  right[n][j], "to bulka lewa")

                if sequence_in_progress != 2: # jak poprzednio
                    struct[n].append({'type':'bulge', 'start_l': i, 'length':1 })
                    sequence_in_progress = 2
                elif sequence_in_progress == 2:
                    struct[n][-1]["length"] += 1
                i += 1 # jak poprzednio
                # print("calutka lista:", struct)
                # print("dlugosc struktury", structure_length)
            else:
                break
    return struct


linear_1st_structures_list = find_slb(left_hairpin_side, right_hairpin_side, structures)
print(linear_1st_structures_list)


# DO TEJ PORY WSZYSTKO ŁADNIE - TRZEBA TYLKO PRZECZYŚCIĆ Z PRINTÓW

###################################################################################

# normalizacja numerów indeksów aby zgadzały się z indeksami wejściowej sekwencji

for k in range(len(linear_1st_structures_list)):
    for structure in linear_1st_structures_list[k]:
        if "start_l" in structure.keys(): # warunek dodany, gdyż spinka ma tylko 1 początek i koniec a bulge mają długość tylko z 1 strony
            if structure['type'] == 'loop': # długość spinki jest liczona "podwójnie" dlatego konieczna jest pętla do poprawy tego
                if structure['length'] % 2 == 0:
                    structure['end_l'] = structure['start_l'] + (structure['length']/2) # jeśli pętla symetryczna koniec to długośc/2
                elif structure['length'] % 2 == 1:
                    structure['end_l'] = structure['start_l'] + ((structure['length']-1)/2) # jeśli pętla niesymetryczna koniec to (długość-1)/2

                if structure['assymetry'] == 'left': # branie pod uwagę asymetryczności pętli - jeden z końców jest dłuższy niż drugi
                    structure['end_l'] += 1 # ASYMETRIA MOŻE BYĆ > 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            else:
                structure['end_l'] = structure['start_l'] + structure['length']
                # dopisanie do kazdego słownika (spinki/pętli/pnia/bulge) indeks zakończenia struktury po lewej stronie
                # dodając wartości początku i długości

            structure['start_l'] = len(left_hairpin_side[k]) - structure['start_l']
            structure['end_l'] = len(left_hairpin_side[k]) - structure['end_l']
            # normalizacja indeksów sekwencji
            # sekwencje lewe były indeksowane na odwróconym stringu, dlatego od jego dł odejmuje się indeks

        if "start_r" in structure.keys(): # analogiczna pętla dla części prawych
            if structure['type'] == 'loop':
                if structure['length'] % 2 == 0:
                    structure['end_r'] = structure['start_r'] + (structure['length']/2)
                elif structure['length'] % 2 == 1:
                    structure['end_r'] = structure['start_r'] + ((structure['length']-1)/2)

                if structure['assymetry'] == 'right':
                    structure['end_r'] += 1
            else:
                structure['end_r'] = structure['start_r'] + structure['length']


            structure['start_r'] = structure['start_r'] + linear_1st_structures_list[k][0]["end"]
            structure['end_r'] = structure['end_r'] + linear_1st_structures_list[k][0]["end"]
            # normalizacja indeksów sekwencji
            # sekwencje prawe zaczynały się w raz z końcem danej spinki, do tej wartości trzeba dodac indeks



print("Po znormalizowaniu indeksów;", linear_1st_structures_list)

#############################################

# zamiana odpowiednich znaków w sekwencji na nazwy struktur : A0, A1 itp
linear_1st_structures_names = []
for p in range(len(linear_1st_structures_list)):
    start_place = linear_1st_structures_list[p][-1]["end_l"]
    end_place = linear_1st_structures_list[p][-1]["end_r"]
    linear_1st_structures_names.append({"name":"A"+str(p), "start":start_place, "end":end_place})
# stworzenie nowej listy w której każda struktura liniowa otrzyma słownik z kluczami: imieniem, początkim i końcem
print(linear_1st_structures_names)

print(dot_bracket_seq)

dot_bracket_seq_2nd = ""

for p in range(len(linear_1st_structures_names)):
    if p == 0:
        dot_bracket_seq_2nd += dot_bracket_seq[:linear_1st_structures_names[p]["start"]] + linear_1st_structures_names[p]["name"] + dot_bracket_seq[linear_1st_structures_names[p]["end"]:linear_1st_structures_names[p+1]["start"]]

    if p != 0 and p != (len(linear_1st_structures_names)-1):
        dot_bracket_seq_2nd += linear_1st_structures_names[p]["name"] + dot_bracket_seq[linear_1st_structures_names[p]["end"]:linear_1st_structures_names[p+1]["start"]]

    if p == (len(linear_1st_structures_names)-1):
        dot_bracket_seq_2nd += linear_1st_structures_names[p]["name"] + dot_bracket_seq[linear_1st_structures_names[p]["end"]:]

# zapisanie nowej wersji oryginalnej sekwencji, gdzie struktury liniowe zostaną zastąpione ich nazwami
# początek sekwencji do startu 1szej (indeks 0) struktury liniowej, nazwa struktury, zakres od jej końca do kolejnej,
# nazwa kolejnej struktury, zakres od jej końca do kolejnej (razy x)
# nazwa ostatniej struktury, zakres od jej końca do końca sekwencji

print("Nowa wersja oryginalnej sekwencji:", dot_bracket_seq_2nd)

###########################################################

# odnajdywanie skrzyżowań - kropki między strukturami liniowymi
crossroad_pat = "\.*" # wyrażenie regularne odnajdujące skrzyżowanie - kropki pomiędzy strukturami liniowymi
crossroad_pattern = re.compile(crossroad_pat) # kompilacja wyrażenia regularnego

#for linear_structure in linear_1st_structures_names:





#########################################################
# CO JEŚLI ASYMETRIA PĘTLI JEST > 1 ??
# wszystkie nowe struktury połączone tylko kropkami zapisujemy jako skrzyzowanie wraz z jego długością (wszystkie kropki po prawej i lewej aż do najbliższego nawiasu
# powtórzeni funkcji find_slb traktując (skrzyżowanie) jako spinkę do włosów
# inne skrzyżowania znalezione w tym czasie traktowane sa jako kolejna struktura liniowa B0 B1 itd
# znowu rozbudowa algorytmu jak wcześniej aż do momentu kiedy będą tylko ... na końcach


# usuwanie już na początku znaków pseudowęzłów
# odszukiwanie indeksów pseudowęzłów i przypisywanie ich do konkretnych struktur

# deal with igraph




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

# listka = [0,1]
# for p in range(len(listka)):
#     print(p)
# print(len(listka)-1)