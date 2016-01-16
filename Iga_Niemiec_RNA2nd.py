import re
in_file = open("C:/Users/Igusia/Documents/example.txt")
in_file2 = in_file.read().replace("\n", "")
input_file = str(in_file2)
#print("Plik wejsciowy: \n", input_file)
#print()
# importowanie pliku wejściowego, usunięcie newline i zamienienie na string

bases = ["A", "C", "U", "T", "N"]
signs = [".", "(", ")", "[", "]", ">", "<", "{", "}"]
# listy zasad azotowych (po nich następuje właściwa sekwencja dot-bracket, oraz znaków w poszukiewanym formacie

def find_seq(input):
    sequence = ""
    for i in range(len(input)):
        if input[i] in bases and input[i+1] in signs and input[i+2] in signs and input[i+3] in signs:
            sequence = (input[i+1:])
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
hairpin_pat = "\(\.+\)" # wyrazenie regularne znajdujace spinki (dowolna ilość "." zakończona "(" z lewej i ")" z prawej strony)
hairpin_pattern = re.compile(hairpin_pat) # kompilacja wyrażenia regularnego

hairpins = []
for m in hairpin_pattern.finditer(dot_bracket_seq): # pętla dla każdego patternu wyszukanego w sekwencji
    hairpins.append({'start': m.start() - 1, 'end': m.end() - 1, 'length': m.end() - m.start()})
print("Znalezione struktury spinek:\n", hairpins)
# stworzenie listy słowników z których każdy będzie zawierał inf o kolejnych spinkach: rozpoczęcie, zakończenie i długość

structures = []
for n in range(len(hairpins)):
    structures.append([])
    structures[n].append(hairpins[n])
print(structures) # lista z poszczególnymi strukturami liniowymi.
# każda str jest zapisana w osobnej liście skł się ze słowników z informacjami
# o jej składnikach tj spinkach/pniach/pętlach/bulge, miejscach ich rozpoczecia i zakończenia w str i długości


# wyszukiwanie pnia/bulge/pętli wew

left_hairpin_side = []
right_hairpin_side = []
for n in range(0, len(hairpins)):
    print("Spinka ",n, ":", hairpins[n])
    left_side = dot_bracket_seq[0: hairpins[n]['start']]
    # print("Lewa strona od spinki", left_side)
    left_hairpin_side.append(left_side) # lista sekwencji po lewej stronie każdej ze znalezionch spinek
    # print("Lista sekwencji po lewej stronie spinki:", left_hairpin_side)
    right_side = dot_bracket_seq[hairpins[n]['end']: -1]
    # print("Prawa strona od spinki", right_side)
    right_hairpin_side.append(right_side) # analogiczna lista z sekwencjami na prawo od spinek
    # print("Lista sekwencji po prawej stronie spinki:", right_hairpin_side)

# dla każdej znalezionej spinki wypisanie fragmentu sekwencji po jej lewej oraz prawej stronie
# bez uwzględnienia znaków wchodzących w skłąd spinki
# stworzenie list sekwencji "prawostronnych" oraz "lewostronnnych"

print("Lista sekwencji po lewej stronie spinek:", left_hairpin_side)
print("Lista sekwencji po prawej stronie spinek:", right_hairpin_side)


# def find_slb(left, right): # funkcja odnajdująca pnie/pętle/bulge wew(stems/bulges/loops) połączone liniowo ze spinkami
#
#     for n in range(len(left)): # zakres odp ilości "lewych części spinek" a więc i znalezionych spinek
#         print("Dla ", n, "spinki część lewostronna ma długość: ", len(left[n]), " a prawostronna ma długość: ", len(right[n]))
#         structure_lenght = len(structures[n]) # zmienna oznaczająca aktualną długość listy odpowiadającej strukturze liniowej do której dopisywane są kolejne el: pnie/pętle itp
#
#         # for i in range(len(left[n])): # zakres odp długości n-tej "lewej sekwencji"
#         #     print(i)
#         #     sequence_in_progress = 0
#         #
#         #     for j in range(len(right[n])): # dla zakresu znaków znajdujących się w mniejszej części (prawej lub lewej) dla danej sekwencji
#         #         print(j)
#         #         if left[n][-(i+1)] == "(" and right[n][j] == ")": #jeśli po obu stronach spinki nawiasy - identyfikacja jako pień
#         #             print("to pien")
#         #             # if sequence_in_progress == 0:
#         #             #     print(left[n][-(i+1)], "oraz",  right[n][i], "to stem")
#         #             #     structures[n].append({"type":"stem", "startl": left[n][-(i+1)], "startr": right[n][i]})
#         #             #     structure_lenght += 1
#         #             # sequence_in_progress = 1
#         #         else:
#         #             print("to NIE pien")
#     #                 sequence_in_progress = 0
#     #                 structures[n][structure_lenght]["endl"] = -(i+1)
#     #                 structures[n][structure_lenght]["endr"] = i
    # print(structures)

                # if left[n][-(i+1)] == "." and right[n][i] == ".": #jeśli po obu stronach spinki kropki - identyfikacja jako pętla
                #     print(left[n][-(i+1)], "oraz",  right[n][i], "to loop")
                # else:
               # if (left[n][-(i+1)] == "." and right[n][i] == ")") or (left[n][-(i+1)] == "(" and right[n][i] == "."):
        #             print(left[n][-(i+1)], "oraz",  right[n][i], "to bulge")
        # print("Koniec struktury")


### ..(((..(((...)))))))....)



# trzeba jeszcze:
# bulka ma zmieniać liczenie tak jakby jej nie było
# ze znalezieniem struktury dopisanie nowego słownika i jego cech do których będzie dodawane +1 w każdej iteracji
# nowa struktura za każdy razem kiedy poprzednia jest "przerwana"
# koniec kiedy kolejnymi znakami są .. i nawiasy w 2 stronę


# każdy s/b/l
# lista_lewa = ["(((...(((..(", "((((....((((.", "(((..(((.."]
# lista_prawa = [")))....)))..", "))..", "....))))..)"]
#
# print(find_slb(lista_lewa, lista_prawa))

# print(" A teraz nasza sekw:", find_slb(left_hairpin_side, right_hairpin_side))





# re.findall(str[0:pierwsza_spinka_init])
# re.findall(str[pierwsza_spinka_end:end])
# lista_pni[[[12,15],[18,21]],[[]]]


# lista_list = [[{1:'a', 2:'b'}, {3:'c', 4:'d'}], [{'z':24, 'y':23}, {'x':22, 'w':21}], ["hey", "you"]]
# lista_list[0].append({'f':17})
# print("Lista", lista_list)
# print(lista_list[2][0][1])


# print(len(left_hairpin_side[-1]))

# for i in range(0,100):
#     if i in range(0,5) and i in range(0,3):
#         print(i)

print("TESTTESTTEST")
kkk = [1,2,3,4,5]

def dodawanie(l):
    for i in l:
        i = 1
    return l

lll = dodawanie(kkk)
print(lll)