#!/bin/bash 

awk '
BEGIN{
    id = 1
}
{
    for(j=1; j<=NF; j+=1)
    {
        data[id] = $j
        id += 1
    }
}

function prtShape(flagNum, flagName, cellNum, i, data)
{
    if(data[i] == flagNum)
    {   
        printf "%s->[", flagName
        for(j=i+1; j<=i+cellNum; j+=1)
        {
            printf "%d,", data[j]
        }
        printf "]\n"
        i = i + cellNum + 1
    }

    return i
}

END{
    i=1
    while(i<id)
    {
        i = prtShape(5, "TRI_3", 3, i, data)
        i = prtShape(7, "QUAD_4", 4, i, data)
        i = prtShape(10, "TETRA_4", 4, i, data)
        i = prtShape(12, "PYRA_5", 5, i, data)
        i = prtShape(14, "PENTA_6", 6, i, data)
        i = prtShape(17, "HEXA_8", 8, i, data)
    }
}
-f ' $1