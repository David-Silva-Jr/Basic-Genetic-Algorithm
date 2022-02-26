#include <iostream>
#include <exception>
#include "StringArithmetic.h"
#include <vector>
#include <algorithm>
#include <cmath>

struct Gene{
    char name;
    int val;

    Gene(char _name, int _val){
        name = _name;
        val = _val;
    }

    Gene(){
        name = '_';
        val = 0;
    }
};

struct Chromosome{
    int size;
    std::vector<Gene> genes;

    Chromosome(int _size){
        size = _size;
    }

    Chromosome(){
        size = 0;
    }

    Gene* operator[](int _i){
        return &genes[_i];
    }

    Gene* operator[](char _c){
        for(int i = 0; i < size; i++){
            if(genes[i].name == _c){
                return &genes[i];
            }
        }
        return nullptr;
    }
};

struct Solution{
    size_t fitness;
    Chromosome data;

    Solution(){
        fitness = 0;
    }

    Solution(Chromosome& _data){
        data = _data;
        fitness = 0;
    }

    bool operator<(const Solution& _b){
        return fitness < _b.fitness;
    }

    bool operator>(const Solution& _b){
        return fitness > _b.fitness;
    }
};

/**
 * @brief Swap gene names in an expression with the associated values
 * @param _chrom The chromosome that stores the values to be used
 * @param _cs The equation with gene names as variables
 * @return String where variable names are replaced by values
 */
std::string SwapGenesWithValues(Chromosome& _chrom, std::string _cs){
    std::string out = "";
    for(int i = 0; i < _cs.length(); i++){
        char c = _cs[i];
        //std::cout << "\nscanning: " << c << "\n";
        if( isalpha(c) ){
            //std::cout << c << " is alpha\n";
            if(_chrom[c] != nullptr){
                //std::cout << c << " is a gene name\n";
                out += "[" + std::to_string( _chrom[c]->val) + "]";
                //std::cout << c << " value: " << std::to_string( _chrom[c]->val) << "\n";
            }
            else{
                std::cout << "Gene at " << c << " is null\n";
            }
        }
        else{
            out += c;
        }
    }
    
    // std::cout << _cs << "\n";
    // std::cout << out << "\n";
    return out;
}

/**
 * @brief Get how many constraints have been satisfied
 * @param _sol Solution to be tested for fitness
 * @param _constraints Vector of constraints as strings
 * @return Number of constraints satisfied by _sol
 */
size_t GetFitness(Solution& _sol, std::vector<std::string>& _constraints){
    int counter = 0;
    int numLeft = _constraints.size(); // Allows me to weigh the fitness by constraints satisfied
    std::string curCs;

    for(std::string cs : _constraints){
        //std::cout << cs << "\n";
        //curCs = SwapGenesWithValues(_sol.data, cs);
        curCs = InfixToPostfix(cs);
        //std::cout << curCs << "\n";
        //curCs = InfixToPostfix(curCs);
        curCs = SwapGenesWithValues(_sol.data, curCs);
        //std::cout << curCs << "\n";
        counter += EvaluatePostEquation(curCs);
        if(counter == 0){
            numLeft--;
        }
        //std::cout << std::to_string(counter) << "\n";
    }

    return counter;
}

/**
 * @brief Get the "error" or "cost" of a solution, measured in total difference between the left and right-hand sides of all constraint equations, multiplied by the product of the number of constraints not satisfied and some weight.
 * @param _sol The solution to calculate cost for 
 * @param _constraints Vector of constraints as strings
 * @return Total error
 */
size_t GetError(Solution& _sol, std::vector<std::string>& _constraints){
    size_t error = 0;
    int numLeft = _constraints.size(); // Allows me to weigh the fitness by constraints satisfied
    std::string curCs;

    // Convert constraints into a usable form then return their values
    for(std::string cs : _constraints){
        curCs = InfixToPostfix(cs);
        curCs = SwapGenesWithValues(_sol.data, curCs);
        size_t added = EvaluatePostEquationError(curCs);
        error += added;

        // If constraint is satisfied, deduct 1 from remaining contraints
        if(added == 0){
            numLeft--;
        }
    }

    // Number of constraints not satisfied gets added to the cost
    // This number is weighted heavily to encourage satisfaction 
    // of constraints over just having a low error
    return error+250*numLeft;
}

/**
 * @brief Get total fitness of a vector of solutions
 * @param _sols Vector of solutions
 * @param _constraints Vector of constraints as strings
 * @return Total fitness
 */
size_t TotalFitness(std::vector<Solution>& _sols){
    int count = 0;
    for(int i=0; i < _sols.size(); i++){
        //std::cout << "Fitness of sol " << std::to_string(i) << ": " << std::to_string(_sols[i].fitness) << "\n";
        count += _sols[i].fitness;
    }
    return count;
}

/**
 * @brief Prints information for one solution, such as gene names, values, and fitness
 * @param _sol The solution to print
 */
void PrintSolution(Solution& _sol){

    for(int j = 0; j < _sol.data.size; j++){
        std::cout << _sol.data[j]->name << ": " << std::to_string(_sol.data[j]->val) << " ";
    }
    std::cout << "Error: " << std::to_string(_sol.fitness);
    std::cout << "\n";

}

/**
 * @brief Prints information for a vector of solutions, such as gene names, values, and fitness
 * @param _sol The solutions to print, as a vector
 */
void PrintSolutions(std::vector<Solution>& _sols){
    for(int i=0; i < _sols.size(); i++){
        std::cout << "Solution " << std::to_string(i) << ": ";

        for(int j = 0; j < _sols[j].data.size; j++){
            std::cout << _sols[i].data[j]->name << ": " << std::to_string(_sols[i].data[j]->val) << " ";
        }
        std::cout << "Fitness: " << std::to_string(_sols[i].fitness);
        std::cout << "\n";
    }
}

/**
 * @brief Clips values so they stay in range, also increments nva, should be called whenever a gene is assigned to
 * @todo Possibly move this into the actual gene struct
 * @param _val Value you want to assign to the gene
 * @param _low Lower bound for _val
 * @param _high Upper bound for _gene
 * @param _nva A reference to the nva counter
 */
void EnforceBounds(int& _val, int _low, int _high, size_t& _nva){
    (_val < _low) ? _val = _low : _val = _val;
    (_val > _high) ? _val = _high : _val = _val;

    // I call this every time a gene value is assigned
    _nva++;
}

/**
 * @brief Take _k solutions and return the best one
 * @details Randomly select _k solutions from a vector, with replacement, and return the one with the lowest error
 *          Should be used when fitness has been calculated from error
 * @param _sols Vector of solutions to pick from
 * @param _k Number of solutions to pick
 * @return The solution with the lowest error of those picked
 */
Solution SelectTourneyError(std::vector<Solution>& _sols, int _k){
    std::vector<Solution> selected;

    for(int i = 0; i < _k; i++){
        // Take from anywhere in the vector
        selected.push_back(_sols[rand()%_sols.size()]);
    }
    std::sort(selected.begin(), selected.end());

    // Return the highest value one
    return selected[0];
}

/**
 * @brief Take _k solutions and return the best one
 * @details Randomly select _k solutions from a vector, with replacement, and return the one with the most constraints satisfied
 *          Should be used when fitness has been calculated from constraints satisfied
 * @param _sols Vector of solutions to pick from
 * @param _k Number of solutions to pick
 * @return The solution with the most satisfied constraints of those picked
 */
Solution SelectTourneySatisfied(std::vector<Solution>& _sols, int _k){
    std::vector<Solution> selected;

    for(int i = 0; i < _k; i++){
        selected.push_back(_sols[rand()%_sols.size()]);
    }
    std::sort(selected.begin(), selected.end());

    return selected[selected.size()-1];
}

/**
 * @brief Take two parent solutions and return a child solution
 * @details The child solution will have genes from both parents, with a chance of mutation
 * @param _p1 Parent solution 1
 * @param _p2 Parent solution 2
 * @param _mChance Chance that a gene will mutate on any given assignment
 * @param _mDelta The size of the mutation. The value of a mutation is in the range [-_mDelta, _mDelta]
 * @param _nva Reference to the nva counter
 * @return The child solution
 */
Solution GetOffspring(Solution& _p1, Solution& _p2, float _mChance, int _mDelta, size_t& _nva){
    // Get number of genes
    int numGenes = _p1.data.size;

    // Create chromosome of proper size
    Chromosome chr(numGenes);

    // float to store mutation rolls
    float mutRoll;

    // Do recombination/mutation
    for(int i = 0; i < numGenes; i++){

        // Get the proper data
        char name = _p1.data[i]->name;
        int valP1 = _p1.data[i]->val;
        int valP2 = _p2.data[i]->val;

        // Get value of random parent gene
        int endVal = ((rand()%2) == 0 ? valP1 : valP2);

        // Does a mutation happen?
        mutRoll = (float)rand()/RAND_MAX;
        if(mutRoll < _mChance){
            // Add a number in [-_mDelta, _mDelta] to endVal
            int delta = rand()%(2*_mDelta + 1) - _mDelta;
            endVal += delta;
        }

        // Make sure endVal is in a valid range
        EnforceBounds(endVal, 1, 100, _nva);

        // Add gene to chromosome
        chr.genes.push_back(Gene(name, endVal));
    }

    // return child solution
    return Solution(chr);
}

/**
 * @brief Replaces current generation with a descendant generation
 * @details This version calculates fitness from constraints satisfied
 *          Parent selection is done by Tourney Selection
 * @param _sols Current generation as vector of solutions
 * @param _constraints Vector of constraints to compute fitness against
 * @param _tourneyK Number of contestants per tourney
 * @param _mChance Chance that a gene will mutate on any given assignment
 * @param _mDelta The size of the mutation. The value of a mutation is in the range [-_mDelta, _mDelta]
 * @param _nva Reference to the nva counter
 */
void RepopulateSatisfied(std::vector<Solution>& _sols, std::vector<std::string> _constraints, int _tourneyK, float _mChance, int _mDelta, size_t& _nva){
    std::vector<Solution> outVec;

    for(int i=0; i < _sols.size(); i++){
        //std::cout << "Beginning iteration " << std::to_string(i) << "\n";     
        // I define these in a later line to make it a little easier to try different techniques   
        Solution s1;
        Solution s2;

        // Do tournament selection
        s1 = SelectTourneySatisfied(_sols, _tourneyK);
        s2 = SelectTourneySatisfied(_sols, _tourneyK);

        //std::cout <<  "selected parents\n";

        //PrintSolutions(s1);

        // Create child solution
        Solution s = GetOffspring(s1, s2, _mChance, _mDelta, _nva);
        //std::cout << "Created\n";

        //PrintSolutions(s);

        // Calculate fitness of solution
        //s.fitness = GetDistance(s, _constraints);

        s.fitness = GetFitness(s, _constraints);
        //std::cout << "Got fitness\n";

        // Add solution to population
        outVec.push_back(s);
        //std::cout << "Solution pushed\n";
    }

    // Set active population
    _sols = outVec;
}

/**
 * @brief Replaces current generation with a descendant generation
 * @details This version calculates fitness from error
 *          Parent selection is done by Tourney Selection
 * @param _sols Current generation as vector of solutions
 * @param _constraints Vector of constraints to compute fitness against
 * @param _tourneyK Number of contestants per tourney
 * @param _mChance Chance that a gene will mutate on any given assignment
 * @param _mDelta The size of the mutation. The value of a mutation is in the range [-_mDelta, _mDelta]
 * @param _nva Reference to the nva counter
 */
void RepopulateError(std::vector<Solution>& _sols, std::vector<std::string> _constraints, int _tourneyK, float _mChance, int _mDelta, size_t& _nva){
    std::vector<Solution> outVec;

    for(int i=0; i < _sols.size(); i++){
        // I define these in a later line to make it a little easier to try different techniques   
        Solution s1;
        Solution s2;

        // Do tournament selection
        s1 = SelectTourneyError(_sols, _tourneyK);
        s2 = SelectTourneyError(_sols, _tourneyK);

        // Create child solution
        Solution s = GetOffspring(s1, s2, _mChance, _mDelta, _nva);

        // Calculate fitness of solution
        s.fitness = GetError(s, _constraints);

        // Add solution to population
        outVec.push_back(s);
    }

    // Set active population
    _sols = outVec;
}

/**
 * @brief Create new population from single individual
 * 
 * @param _sols Vector of solutions to be replaced
 * @param _constraints Set of constraints to measure fitness against
 * @param _mChance Chance of mutation during replication of a given gene
 * @param _mDelta Size of mutation
 * @param _nva Variable assignment counter
 */
void RepopulateFromProgenitor(std::vector<Solution>& _sols, std::vector<std::string> _constraints, float _mChance, int _mDelta, size_t& _nva){
    std::vector<Solution> outVec;
    std::sort(_sols.begin(), _sols.end());

    for(int i=0; i < _sols.size(); i++){
        // Create child solutions
        Solution s = GetOffspring(_sols[0], _sols[0], _mChance, _mDelta, _nva);

        // Calculate fitness of solution
        s.fitness = GetError(s, _constraints);

        // Add solution to population
        outVec.push_back(s);
    }

    // Set active population
    _sols = outVec;
}

/**
 * @brief Finds the solution with the lowest error
 * @param _sols Vector of solutions
 * @return The solution with the lowest error
 */
Solution& BestSolutionError(std::vector<Solution> _sols){
    // Get ready to loop
    int size = _sols.size();
    int minVal = INT_MAX;
    int minIndex = 0;

    for(int i = 0; i < size; i++){
        // Alias _sols[i] to something easier to read
        Solution& current = _sols[i];
        if(current.fitness < minVal){
            minVal = current.fitness;
            minIndex = i;
        }
    }

    return _sols[minIndex];
}

/**
 * @brief Finds the solution with the most constraints satisfied
 * @param _sols Vector of solutions
 * @return The solution which satisfied most constraints
 */
Solution& BestSolutionSatisfied(std::vector<Solution> _sols){
    // Get ready to loop
    int size = _sols.size();
    int most = -1;
    int mostIndex = 0;
    for(int i = 0; i < size; i++){
        // Alias to make this easier to read
        Solution& current = _sols[i];
        if(current.fitness > most){
            most = current.fitness;
            mostIndex = i;
        }
    }

    return _sols[mostIndex];
}

/**
 * @brief Print a list of constraints, whether or not they're satisfied, and how far they are from satisfaction based on their error value
 * @param _sol Solution to be tested
 * @param _constraints List of contraints to test against
 */
void ReportError(Solution _sol, std::vector<std::string> _constraints){
    std::string curCs;

    for(std::string cs : _constraints){
        // Convert and compute constraint
        curCs = InfixToPostfix(cs);
        curCs = SwapGenesWithValues(_sol.data, curCs);
        int error = EvaluatePostEquationError(curCs);

        // Checkbox is [x] if satisfied and [ ] if not
        std::cout << "      [";
        if(error == 0){
            std::cout << 'x';
        }
        else{
            std::cout << " ";
        }
        std::cout<< "]";

        // Print unaltered constraint
        std::cout << cs << "\n";

        // Print constraint with proper values and the error
        std::cout << "         " << SwapGenesWithValues(_sol.data, cs) << ", error: ";
        std::cout << std::to_string(error) << "\n";
    }
}

/**
 * @brief Print a list of constraints, whether or not they're satisfied
 * @param _sol Solution to be tested
 * @param _constraints List of contraints to test against
 */
void ReportSatisfied(Solution _sol, std::vector<std::string> _constraints){
    std::string curCs;

    for(std::string cs : _constraints){
        curCs = InfixToPostfix(cs);
        curCs = SwapGenesWithValues(_sol.data, curCs);
        bool satisfied = EvaluatePostEquation(curCs);

        std::cout << "      [";
        if(satisfied){
            std::cout << 'x';
        }
        else{
            std::cout << " ";
        }
        std::cout<< "]";
        std::cout << cs << "\n";
        std::cout << "         " << SwapGenesWithValues(_sol.data, cs);
    }
}

/**
 * @brief Compute rate of mutation at a given time
 * @param _begin Rate of mutation at time 0
 * @param _end Target rate of mutation
 * @param _rate The rate at which the function tends towards _end
 * @param _t The current time, i.e. the current generation count
 * @return Float representing the current mutation rate
 */
float GetMutationRateAtTime(float _begin, float _end, float _rate, int _t ){
    return _end + (_begin - _end) * std::exp(-1.0*_rate*_t);
}

/**
 * @brief Performs genetic computing to solve a set of constraints. Fitness based on solutions satisfied.
 * 
 * @param _geneCount Number of genes per chromosome
 * @param _populationSize Number of solutions per generation
 * @param _mutationDelta How large a mutation can be
 * @param _tourneyK Number of contestants for Tourney Selection
 * @param _gene0 Char representing the name of the first gene (If this is 'A,' then subsequent genes will be 'B,' 'C,' etc.)
 * @param _baseMutationChance At time 0, how likely is it that a gene will mutate on replication
 * @param _endMutationChance Probability that a gene will mutate on replication at some point in the future
 * @param _coolRate Rate at which chance of mutation goes from base to end
 * @param _solutions Vector of solutions, can be empty as the function initializes them
 * @param _constraints List of constraints as a vector of strings
 * @param _nva Counter for number of variable assignments
 */
void DoGenerationsSatisfied(int _geneCount, int _populationSize, int _mutationDelta, int _tourneyK, char _gene0, float _baseMutationChance, float _endMutationChance, float _coolRate, std::vector<Solution>& _solutions, std::vector<std::string> _constraints, size_t& _nva){
    // Initial generation with constraints satisfied fitness
    for(int i=0; i < _populationSize; i++){
        Chromosome data(_geneCount);
        for(int j = 0; j < _geneCount; j++){
            data.genes.push_back(Gene(_gene0 + j, (rand()%100)+1));
        }
        _solutions.push_back(Solution(data));
        _solutions[i].fitness = GetFitness(_solutions[i], _constraints);
    }

    Solution best = BestSolutionSatisfied(_solutions);
    float coolFactor;
    int i = 0;
    while(best.fitness < _constraints.size() && i < 3000){
        // Compute current cooling factor
        coolFactor = GetMutationRateAtTime(_baseMutationChance, _endMutationChance, _coolRate, i);
        //std::cout << "Cooling factor: " << std::to_string(coolFactor) << "\n";


        // Replace population
        RepopulateSatisfied(_solutions, _constraints, _tourneyK, coolFactor, _mutationDelta, _nva);
        
        // Get best after doing repopulation
        best = BestSolutionSatisfied(_solutions);

        // Report every now and then
        if(i%100 == 0){
            std::cout << "\n----------------------------GENERATION " << std::to_string(i) << " (" << std::to_string(_constraints.size()) << " constraints) " << "----------------------------------\n";
            std::cout << "      Mutation chance: " << std::to_string(100*coolFactor) << "%\n";
            std::cout << "      Average error: " << std::to_string((float)TotalFitness(_solutions)/_solutions.size()) << "\n";
            std::cout << "      Best error: ";
            PrintSolution(best);
            ReportSatisfied(best, _constraints);
            //PrintSolutions(solutions);
        }

        i++;
    }

    if(i >= 3000){
        std::cout << "No solution found within 3000 generations.\n";
        return;
    }

    // Report solutions
    std::sort(_solutions.begin(), _solutions.end());
    PrintSolutions(_solutions);
    PrintSolution(best);
    ReportSatisfied(best, _constraints);
    std::cout << "Solved in " << std::to_string(i) << " generations.\n";
}

/**
 * @brief Performs genetic computing to solve a set of constraints. Fitness based on error.
 * 
 * @param _geneCount Number of genes per chromosome
 * @param _populationSize Number of solutions per generation
 * @param _mutationDelta How large a mutation can be
 * @param _tourneyK Number of contestants for Tourney Selection
 * @param _gene0 Char representing the name of the first gene (If this is 'A,' then subsequent genes will be 'B,' 'C,' etc.)
 * @param _baseMutationChance At time 0, how likely is it that a gene will mutate on replication
 * @param _endMutationChance Probability that a gene will mutate on replication at some point in the future
 * @param _coolRate Rate at which chance of mutation goes from base to end
 * @param _solutions Vector of solutions, can be empty as the function initializes them
 * @param _constraints List of constraints as a vector of strings
 * @param _nva Counter for number of variable assignments
 */
void DoGenerationsError(int _geneCount, int _populationSize, int _mutationDelta, int _tourneyK, char _gene0, float _baseMutationChance, float _endMutationChance, float _coolRate, std::vector<Solution>& _solutions, std::vector<std::string> _constraints, size_t& _nva){
    // Initial generation with error fitness
    for(int i=0; i < _populationSize; i++){
        Chromosome data(_geneCount);
        for(int j = 0; j < _geneCount; j++){
            data.genes.push_back(Gene(_gene0 + j, (rand()%100)+1));
            _nva++;
        }
        _solutions.push_back(Solution(data));
        _solutions[i].fitness = GetError(_solutions[i], _constraints);
    }

    // Store the best solution
    Solution best = BestSolutionError(_solutions);
    float mutationFactor;
    int i = 0;
    while(best.fitness > 0 && i < 3000){
        // Compute current mutation rate
        mutationFactor = GetMutationRateAtTime(_baseMutationChance, _endMutationChance, _coolRate, i);

        // Replace population
        RepopulateError(_solutions, _constraints, _tourneyK, mutationFactor, _mutationDelta, _nva);
        
        // Get best after doing repopulation
        best = BestSolutionError(_solutions);

        // Report every now and then
        if(i%100 == 0){
            std::cout << "\n----------------------------GENERATION " << std::to_string(i) << " (" << std::to_string(_constraints.size()) << " constraints) " << "----------------------------------\n";
            std::cout << "      Mutation chance: " << std::to_string(100*mutationFactor) << "%\n";
            std::cout << "      Average error: " << std::to_string((float)TotalFitness(_solutions)/_solutions.size()) << "\n";
            std::cout << "      Best error: ";
            PrintSolution(best);
            ReportError(best, _constraints);
        }

        i++;
    }

    if(i >= 3000){
        std::cout << "No solution found within 3000 generations.\n";
        return;
    }

    // Report solutions
    std::sort(_solutions.begin(), _solutions.end());
    PrintSolutions(_solutions);
    PrintSolution(best);
    ReportError(best, _constraints);
    std::cout << "Solved in " << std::to_string(i) << " generations.\n";
}

/**
 * @brief Performs genetic computing to solve a set of constraints, working backwards from the highest to the lowest levels.
 * 
 * @param _geneCount Number of genes per chromosome
 * @param _populationSize Number of solutions per generation
 * @param _mutationDelta How large a mutation can be
 * @param _tourneyK Number of contestants for Tourney Selection
 * @param _gene0 Char representing the name of the first gene (If this is 'A,' then subsequent genes will be 'B,' 'C,' etc.)
 * @param _baseMutationChance At time 0, how likely is it that a gene will mutate on replication
 * @param _endMutationChance Probability that a gene will mutate on replication at some point in the future
 * @param _coolRate Rate at which chance of mutation goes from base to end
 * @param _solutions Vector of solutions, can be empty as the function initializes them
 * @param _constraints List of constraints as a vector of strings
 * @param _nva Counter for number of variable assignments
 */
void DoGenerationsBackwards(int _geneCount, int _populationSize, int _mutationDelta, int _tourneyK, char _gene0, float _baseMutationChance, float _endMutationChance, float _coolRate, std::vector<Solution>& _solutions, std::vector<std::string> _constraints, size_t& _nva){
    // Small constraint vector to store constraints currently being worked on
    std::vector<std::string> set = {_constraints[_constraints.size()-1]};

    // Initial generation with error fitness
    for(int i=0; i < _populationSize; i++){
        Chromosome data(_geneCount);
        for(int j = 0; j < _geneCount; j++){
            data.genes.push_back(Gene(_gene0 + j, (rand()%100)+1));
            _nva++;
        }
        _solutions.push_back(Solution(data));
        _solutions[i].fitness = GetError(_solutions[i], set);
    }

    // Store the best solution
    Solution best = BestSolutionError(_solutions);
    float mutationFactor;
    int i = 0;
    while(i < 100000){
        // Compute current mutation rate
        mutationFactor = GetMutationRateAtTime(_baseMutationChance, _endMutationChance, _coolRate, i);

        // Replace population
        RepopulateError(_solutions, set, _tourneyK, mutationFactor, _mutationDelta, _nva);
        
        // Get best after doing repopulation
        best = BestSolutionError(_solutions);

        // Report every now and then
        if(i%100 == 0){
            std::cout << "\n----------------------------GENERATION " << std::to_string(i) << " (" << std::to_string(_constraints.size()) << " constraints) " << "----------------------------------\n";
            std::cout << "      Mutation chance: " << std::to_string(100*mutationFactor) << "%\n";
            std::cout << "      Average error: " << std::to_string((float)TotalFitness(_solutions)/_solutions.size()) << "\n";
            std::cout << "      Best error: ";
            PrintSolution(best);
            ReportError(best, set);
        }

        i++;

        if(best.fitness == 0){
            if(set.size() <= _constraints.size()){
                std::cout << "\n\nFound solution to a new set of constraints. Restarting evolution process\n\n";
                PrintSolution(best);
                ReportError(best, set);
                set.push_back(_constraints[_constraints.size() - set.size() - 1]);
                RepopulateFromProgenitor(_solutions, set, _baseMutationChance, _mutationDelta, _nva);
                i = 0;
            }
            else{
                break;
            }
        }
    }

    if(i >= 100000){
        std::cout << "No solution found within 100000 generations.\n";
        return;
    }

    // Report solutions
    std::sort(_solutions.begin(), _solutions.end());
    PrintSolutions(_solutions);
    PrintSolution(best);
    ReportError(best, _constraints);
    std::cout << "Solved in " << std::to_string(i) << " generations.\n";
}


int main(int argc, char** argv){
    std::vector<Solution> solutions;    // The population of solutions

    int geneCount = 15;                 // How many genes per chromosome?
    int populationSize = 50;            // How many solutions per population?
    int mutationDelta = 1;              // How much does a mutation change the value?
    int tourneyK = 3;                   // How many solutions are selected for tournament selection?
    int seed = 628;                     // for reprodicibility
    int part;                           // Which set of constraints are you using?

    char gene0 = 'A';                   // I didn't want to type 15 different letters

    float baseMutationChance = .95;     // Mutation chance at generation 0
    float endMutationChance = .05;      // Mutation chance much later
    float coolRate = .05;               // Rate of cooling, made negative in function, so must be positive

    // How many times variables have been assigned
    size_t nva = 0;

    std::vector<std::string> constraintsA = {
        // Part A constraints
        "A=B+C+E+F",
        "D=E+F+21",
        "D^2=E*E*A+694",
        "E+F<A",
    };

    std::vector<std::string> constraintsB = {
        // Part A constraints
        "A=B+C+E+F",
        "D=E+F+21",
        "D^2=E*E*A+694",
        "E+F<A",

        // Part B constraints
        "(H*J)+(E*16)=(G+I)^2-48",
        "A-C=(H-F)^2+4",
        "4*J=G^2+7",
    };

    std::vector<std::string> constraintsC = {
        // Part A constraints
        "A=B+C+E+F",
        "D=E+F+21",
        "D^2=E*E*A+694",
        "E+F<A",

        // Part B constraints
        "(H*J)+(E*16)=(G+I)^2-48",
        "A-C=(H-F)^2+4",
        "4*J=G^2+7",

        "O*F*A<M^2",
        "K<14",
        "L^3<M^2",
        "L*O-1<A^2",

        // Part C constraints
        "2*M=K^2-25",
        "(N-O)^2=(J-F)*O*2",
        "N^2=M*J+100",
        "(L+N)^2+1875^2=G*(B+F)*(K+M+N+30)",
        "L*O=(A^2)*(K-G)",
        "L^3=M^2-(O*F*A)",
    };

    // Take command line input
    populationSize = std::stoi(argv[1]);
    mutationDelta = std::stoi(argv[2]);
    tourneyK = std::stoi(argv[3]);
    coolRate = (float)1/std::stoi(argv[4]);
    seed = std::stoi(argv[5]);
    part = std::stoi(argv[6]);

    // Set random seed
    srand(seed);

    // Select constraints
    std::vector<std::string>* selectedConstraints;
    if(part==1){
        selectedConstraints = &constraintsA;
    }
    else if(part==2){
        selectedConstraints = &constraintsB;
    }
    else{
        selectedConstraints = &constraintsC;
    }

    std::cout << "Doing evolution with population: " << std::to_string(populationSize) << ", mutation size: " << std::to_string(mutationDelta) << ", tourney K: " << std::to_string(tourneyK) << ", slowing rate: " << std::to_string(coolRate) << ", seed: " << std::to_string(seed) << ", and constraint set " << part << "\n\n";

    //DoGenerationsSatisfied(geneCount, populationSize, mutationDelta, tourneyK, gene0, baseMutationChance, endMutationChance, coolRate, solutions, constraints, nva);
    //DoGenerationsError(geneCount, populationSize, mutationDelta, tourneyK, gene0, baseMutationChance, endMutationChance, coolRate, solutions, constraintsC, nva);
    DoGenerationsError(geneCount, populationSize, mutationDelta, tourneyK, gene0, baseMutationChance, endMutationChance, coolRate, solutions, *selectedConstraints, nva);

    // Print nva
    std::cout << "\nnva: " << std::to_string(nva) << "\n";

    
}