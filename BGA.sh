g++ "Task2.cpp" -o Task2

echo "Hello and welcome to my version of Task2."
echo 'Are you checking part A, B, or C?'
echo "Warning: C will run for 3000 generations and then self-terminate"
echo "(Enter '1', '2', or '3' with no quotes. If other input is used you must specify parameters.)"
read -p "Selection: " selected

declare popSize
declare alpha
declare k
declare rate
declare seed
declare part

if [ "$selected" = "1" ]
then
    let popSize="25"
    let alpha="5"
    let k="3"
    let rate="20"
    let seed="628"
    let part="1"
elif [ "$selected" = "2" ] 
then
    let popSize="50"
    let alpha="1"
    let k="2"
    let rate="1000"
    let seed="628"
    let part="2"
elif [ "$selected" = "3" ] 
then
    let popSize="50"
    let alpha="1"
    let k="3"
    let rate="1000"
    let seed="628"
    let part="3"
else
    read -p 'Population size: ' popSize
    read -p 'Size of mutations: ' alpha
    read -p 'Tournament K: ' k
    read -p 'Mutation slowing rate (positive integer; the rate will be the negative inverse of your input): ' rate
    read -p 'Seed: ' seed
    read -p 'Constraints from part (as number):' part
fi

./Task2 $popSize $alpha $k $rate $seed $part
