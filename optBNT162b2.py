# Beispielcode für Aufgabe09
# Frank Zimmermann
# 27.04.2021

import csv
import logging
import sys
import time

start = time.time()
# Logger konfigurieren ========================================================
handler1 = logging.StreamHandler(sys.stdout)
handler2 = logging.FileHandler(filename="optimzation.log")
logfmt1 = logging.Formatter(
    "{asctime:24s} {levelname:7s} {filename:15s} {funcName:20s} [{thread:15d}]: {message}",
    style="{")
logfmt2 = logging.Formatter(
    "{asctime} {levelname} {filename} {funcName} [{thread}]: {message}", style="{")
handler1.setFormatter(logfmt1)
handler2.setFormatter(logfmt2)
logger = logging.getLogger()
logger.addHandler(handler1)
logger.addHandler(handler2)
logger.setLevel(logging.INFO)


# Einlesen der csv ============================================================
def read_csv(filename):
    records = []
    with open(filename, 'rt') as fp:
        reader = csv.reader(fp, delimiter=',')
        for ii, row in enumerate(reader):
            if ii > 0:  # Titelzeile ignorieren
                records.append(row)

    return records


# Beispiel--Algorithmus:
# Wenn die 3. Base eines Codons G oder C ist: mache garnichts
# falls nicht: Versuche G als dritte Base:
#              Wenn die Aminosäure gleich bleibt:
#                   fertig.
#              Falls nicht:
#              Versuche C als dritte Base:
#              Wenn die Aminosäure gleich bleibt:
#                   fertig.
# Prüfen, ob Übereinstimmung mit BNT162b2 und ggfs. zähler inkrementieren
def optimization(virvac, amino):
    m = 0
    for element in virvac:
        vir = element[1]
        vac = element[2]
        logger.debug(f'{vir} v {vac}, amino: {amino[vir]} == {amino[vac]}.')

        our_first = change_first(vir, amino)
        our_last = change_last(vir, amino)
        m = count_changes(m, vac, our_first, our_last)

    return m

# count changes in codon-table

def count_changes(m, vac, our_first, our_last):
    if vac == our_last and vac == our_first:
        logger.debug('Matched the vaccine!')
        m += 2
    elif vac == our_last or vac == our_first:
        logger.debug('Matched the vaccine!')
        m += 1
    else:
        logger.debug('No Match.')
    return m

#change last basetriplett

def change_last(vir, amino):
    our = vir
    if vir[2] == 'G' or vir[2] == 'C':
        logger.debug('codon ended on G or C already, not doing anything')
    else:
        prop = vir[:2] + "G"
        logger.debug(f'Attempting G substution, new candidate {prop}')

        if amino[vir] == amino[prop]:
            logger.debug('amino acid still the same, done!')
            our = prop
        else:
            logger.debug(f'Oops, amino acid changed. Trying C, new candidate {prop}')
            prop = vir[:2] + "C"

            if amino[vir] == amino[prop]:
                logger.debug('Amino acid still the same, done!')
                our = prop
    return our

# change first basetriplett

def change_first(vir, amino):
    our = vir
    if vir[0] == 'G' or vir[0] == 'C':
        logger.debug('codon ended on G or C already, not doing anything')
    else:
        prop = "G" + vir[:2]
        logger.debug(f'Attempting G substution, new candidate {prop}')

        if amino[vir] == amino[prop]:
            logger.debug('amino acid still the same, done!')
            our = prop
        else:
            logger.debug(f'Oops, maino acid changed. Trying C, new candidate {prop}')
            prop = "C" + vir[:2]

            if amino[vir] == amino[prop]:
                logger.debug('Amino acid still the same, done!')
                our = prop
    return our

if __name__ == "__main__":
    codons = read_csv('codon-table-grouped.csv')
    logger.debug(codons)
    c2s = {}        # c2s: key=Codon, value=amino
    for c in codons:
        c2s[c[1]] = c[0]
    logger.debug(c2s)
    logger.debug("Länge von c2s: {}".format(len(c2s)))
    vv = read_csv("side-by-side.csv")
    logger.debug(vv)
    logger.debug("Länge von virvac: {}".format(len(vv)))

    length = len(vv)
    tmax = length
    loops = 20      # statistical smoothing
    sumo = 0        # sum of optimization events (mean)
    for i in range(loops):  # Mittelwert Übereinstimmung
        sumo += optimization(vv, c2s)
    sumo /= loops
    for el in vv:  # Zählen der Aminosäuren (Prolin) Ersetzungen
        if c2s[el[1]] != c2s[el[2]]:
            logger.debug(f'{el[1]} v {el[2]}, amino: {c2s[el[1]]} == {c2s[el[2]]}.')
            tmax -= 1
    print("Maximalwert ohne Prolin--Ersetzungen: {:7.2f} %".format(100 * tmax / length))
    print("Durchschnittlich veränderte Codons  :   {:4d}".format(int(sumo)))
    print("Erreichte Übereinstimmung           : {:7.2f} %".format(100 * sumo / length))
    ende = time.time()
    print('{:5.3f}s'.format(ende-start),end= ' ')