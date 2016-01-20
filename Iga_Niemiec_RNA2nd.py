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

############################

def find_slb(left, right, struct):  # funkcja odnajdująca pnie/pętle/bulge wew(stems/bulges/loops) połączone liniowo ze spinkami
    print("Lewa:", left, "i prawa:", right)



    ending_left_pat = "\.*\)" # wyrażenie regularne odnajdujące koniec lewostronny danej struktury liniowej (tj kropki i nawias w odwrotą stronę)
    ending_left_pattern = re.compile(ending_left_pat) # kompilacja wyrażenia regularnego
    ending_right_pat = "\.*\(" # jw. - koniec prawostronny
    ending_right_pattern = re.compile(ending_right_pat) # kompilacja wyrażenia regularnego


    for n in range(len(left)):  # zakres odp ilości "lewych części spinek" a więc i znalezionych spinek
        rev_left = left[n][::-1] # odwrócenie sekwencji "lewych cześci spinek" aby łatwiejsze było porównywanie sekwencji "prawej" i "lewej"
        print("Odwrotność:", rev_left)
        print("Dla ", n, "spinki część lewostronna ma długość: ", len(left[n]), " a prawostronna ma długość: ", len(right[n]))
        sequence_in_progress = 0  # zmienna oznaczająca strukturę w której pętli obecni znajduje się funkcja
        i = 0 # zmienna oznaczająca znak w sekwencji "lewej"
        j = 0 # zmienna oznaczająca znak w sekwencji "prawej"

        ending_places_left = []
        ending_places_right = []

        for end_l in ending_left_pattern.finditer(rev_left): # pętla dla każdego patternu wyszukanego w sekwencji
            ending_places_left.append(end_l.start())
            print("Znalezione miejsca kończące w lewej części od spinki:\n", ending_places_left)
        for end_r in ending_right_pattern.finditer(right[n]): # pętla dla każdego patternu wyszukanego w sekwencji
            ending_places_right.append(end_r.start())
            print("Znalezione miejsca kończące w prawej części od spinki:\n", ending_places_right)


        if len(ending_places_left) == 0:
            end_place_left = len(rev_left)
        else:
            end_place_left = ending_places_left[0]

        if len(ending_places_right) == 0:
            end_place_right = len(right[n])
        else:
            end_place_right = ending_places_right[0]

        print("KONIEC LEWY", end_place_left, "I PRAWY", end_place_right)


        #while i in range(len(rev_left)) and j in range(len(right[n])): # zakres odp długości " odwróconej lewej" i n-tej "prawej" sekwencji
        while i in range(end_place_left) and j in range(end_place_right):
            print("To jest i:", i, "a to j:", j)
            structure_length = len(struct[n])  # zmienna oznaczająca aktualną długość listy odpowiadającej strukturze liniowej do której dopisywane są kolejne el: pnie/pętle itp

            if rev_left[i] == "(" and right[n][j] == ")":  # jeśli po obu stronach spinki nawiasy - identyfikacja jako pień
                print(rev_left[i], "oraz",  right[n][j], "to pien")

                if sequence_in_progress != 1:
                    struct[n].append({'type':'stem', 'start_l': i, 'start_r': j, 'length':1 })
                    sequence_in_progress = 1 # po wejściu do danej pętli np pnia, do listy dodawany jest nowy słownik zaw inf o strukturze
                elif sequence_in_progress == 1:
                    struct[n][-1]["length"] += 1 # przy kolejnym wejściu w tą samą pętlę dodawane jest jedynie +1 do wartości długości
                i += 1
                j += 1 # przechodzenie symetryczne do kolejnych znaków
                print("calutka lista:", struct)
                print ("Nowe i:", i , "oraz j:", j)
                print("dlugosc struktury", structure_length)

            elif rev_left[i] == "." and right[n][j] == ".": # jeśli po obu stronach spinki kropki - identyfikacja jako pętla
                print(rev_left[i], "oraz",  right[n][j], "to petla")

                if sequence_in_progress != 3:
                    struct[n].append({'type':'loop', 'start_l': i, 'start_r': j, 'length':2 }) # TO TEŻ NIE POWINNO SIĘ TAK DODAWAĆ BO CZASEM JEST NIESYMETRYCZNE
                    sequence_in_progress = 3
                elif sequence_in_progress == 3:
                    struct[n][-1]["length"] += 2 # działanie jak przy poprzedniej pętli


                if rev_left[i+1] == "(" and right[n][j+1] == ".":
                    print(rev_left[i+1], "oraz",  right[n][j+1], "to petla ale niesymetryczna z lewej") # TU TRZEBABY DODAC WARUNEK CO JEŚLI WYJDZIE Z RANGE ALE JEŚLI BEDZIE SIĘ ODCINAĆ
                    # FRAGMENTY JEDNONICIOWE NA KOŃCACH TO NIE TRZEBA
                    j += 1
                elif rev_left[i+1] == "." and right[n][j+1] == ")":
                    print(rev_left[i+1], "oraz",  right[n][j+1], "to petla ale niesymetryczna z prawej")
                    i += 1


                i += 1
                j += 1
                print("calutka lista:", struct)
                print("dlugosc struktury", structure_length)

            elif rev_left[i] == "(" and right[n][j] == ".": # jeśli po obu stronach spinki nawiasy - identyfikacja jako pień
                print(rev_left[i], "oraz",  right[n][j], "to bulka prawa")

                if sequence_in_progress != 2:
                    struct[n].append({'type':'bulge', 'start_r': j, 'length':1 })
                    sequence_in_progress = 2
                elif sequence_in_progress == 2:
                    struct[n][-1]["length"] += 1
                j += 1
                print("calutka lista:", struct)
                print("dlugosc struktury", structure_length)

            elif rev_left[i] == "." and right[n][j] == ")": # jeśli po obu stronach spinki nawiasy - identyfikacja jako pień
                print(rev_left[i], "oraz",  right[n][j], "to bulka lewa")

                if sequence_in_progress != 2:
                    struct[n].append({'type':'bulge', 'start_l': i, 'length':1 })
                    sequence_in_progress = 2
                elif sequence_in_progress == 2:
                    struct[n][-1]["length"] += 1
                i += 1
                print("calutka lista:", struct)
                print("dlugosc struktury", structure_length)

            # teraz będzie kończenie funkcji gdy znajdzie "." i ) po lewej lub ( po prawej

            # elif rev_left[i] == ")" or right[n][j] == "(":
            # #     print("Uff, koniec struktury liniowej - funkcja kończy działanie")
            # #     if struct[n][-1]["type"] == "loop" or struct[n][-1]["type"] == "bulge"
            else:
                break




# trzeba jeszcze:
# koniec kiedy kolejnymi znakami są .. i nawiasy w 2 stronę
# pozmieniać te słowniki: zamienić indeksy z rev_left na left oraz dodać wartości end_l i end_r


# każdy s/b/l
lista_lewa = ["(((...((..(", "((((....((((."]
lista_prawa = [")))....))).((", ")).."]

lista_lewa2 = ['((((((((((.....((((((((....(((((', '((((((((((.....((((((((....(((((((.............))))..)))...)))))).)).(((((((..((((']
lista_prawa2 = ['))))..)))...)))))).)).(((((((..((((((....))))))..)))))))...))))))))))', '))))))..)))))))...))))))))))']

print(find_slb(left_hairpin_side, right_hairpin_side, structures))




#
print()
print("proba dodawania")
lista_list = [[{1:'a', 2:'b'}, {3:'c', 4:'d'}], [{'z':24, 'y':23}, {'x':22, 'w':21}], ["hey", "you"]]
lista_list[0].append({'f':17})
print("Lista", lista_list)
lista_list[1][0]["z"] += 1
print("Lista2", lista_list)
print(lista_list[2][0][1])

def probna_funkcja(lista):
    if lista[1][0]["y"] == 23:
        print("superowo!")
    else:
        print("nie do końca")

probna_funkcja(lista_list)

# print(len(left_hairpin_side[-1]))

# for i in range(0,100):
#     if i in range(0,5) and i in range(0,3):
#         print(i)

k = [1,4,5,6,7,8]
p = []
print(len(p))



