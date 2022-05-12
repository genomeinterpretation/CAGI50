function  [BLUE, GREEN, RED, PURPLE, ORANGE, GREY, DBLUE] = mainColors()
    GREEN = [0.47, 0.67, 0.19];
    RED = [0.85, 0.33, 0.10];
    PURPLE = [0.49, 0.18, 0.56];
    ORANGE = [0.93, 0.69, 0.13];
    GREY = [0.5, 0.5, 0.5];
    BLUE = [0, 0.45, 0.74];
    DBLUE = [39, 50, 112];
    DBLUE = DBLUE/sum(DBLUE);
end