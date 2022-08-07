def getStarColours():
    colours = []
    with open('starColours.txt', mode='r') as file:
        text = file.read().split('\n')
        for line in text:
            r = int(line[13:16])
            g = int(line[17:20])
            b = int(line[21:24])
            colours.append((r, g, b))

    return colours