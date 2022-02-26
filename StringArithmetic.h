#ifndef SARIM
#define SARIM

#include <stack>
#include <string>

// Exponentiation by squaring
// Based on the wikipedia article
// Pretty please don't pass any negative powers in
// @param _base The base of the exponent
// @param _power Power you are raising the base to
int Pow(int _base, int _power){
    if(_power == 0) return 1;
    else if (_power % 2 == 0) return Pow(_base*_base, _power/2);
    else if (_power % 2 == 1) return _base * Pow(_base*_base, (_power-1)/2);
    else return 1;
}


// This function based heavily on the GeeksForGeeks article for postfix evaluation
// https://www.geeksforgeeks.org/stack-set-4-evaluation-postfix-expression/
//
// NOTE: I need to update this or make a separate version that returns true or false
// based on if an input EQUATION is true or not. Also need to make it get the value
// stored in the gene with the same name as the current character and evaluate using
// that.
// @param exp Arithmetic expression in postscript form
// @return value of the expression as an integer
int EvaluatePostfix(std::string exp)
{
    std::stack<int> stack;
 
    // Scan all characters one by one
    for (int i = 0; i < exp.length(); i++)
    {
        char c = exp[i];
        
        // Keep track of where an integer started
        if (c == '['){
            stack.push(c);
        }
        else if (c == ']'){
            // Current sum
            int num = 0;

            // How deep in the stack did I pull the current number from?
            int counter = 0;
            while(stack.top() != '['){
                // Add current number * 10^depth to find actual number value
                num += stack.top() * Pow(10, counter);
                stack.pop();

                // Increment the counter so you know what power to raise the next number to
                counter++;
            }

            // Pop '['
            stack.pop();

            // Push the correct value
            stack.push(num);
        }
        else if (isdigit(c))
            stack.push(c - '0');
 
        // If the scanned character is an operator, pop two
        // elements from stack apply the operator
        else
        {   
            // Remember to pop after getting the value
            int val1 = stack.top();
            stack.pop();

            int val2 = stack.top();
            stack.pop();

            switch (exp[i])
            {
            case '+': stack.push(val2 + val1); break;
            case '-': stack.push(val2 - val1); break;
            case '*': stack.push(val2 * val1); break;
            case '/': stack.push(val2 / val1); break;
            case '^': stack.push(Pow(val2, val1)); break;
            }
        }
    }
    return stack.top();
}

// This function taken from the GeeksForGeeks article for infix to postfix
// https://www.geeksforgeeks.org/stack-set-2-infix-to-postfix/
// @param c Operator to get the precedence of
// @return Precedence of the operator c, or -1 if c is not an operator
int Prec(char c) {
    if(c == '^')
        return 3;
    else if(c == '/' || c=='*')
        return 2;
    else if(c == '+' || c == '-')
        return 1;
    else
        return -1;
}

// This function based heavily on the GeeksForGeeks article for infix to postfix
// https://www.geeksforgeeks.org/stack-set-2-infix-to-postfix/
// The main function to convert infix expression
// to postfix expression
// @param s Arithmetic expression in infix notation to be converted
std::string InfixToPostfix(std::string s) {
 
    std::stack<char> st; //For stack operations, we are using C++ built in stack
    std::string result;
 
    for(int i = 0; i < s.length(); i++) {
        char c = s[i];
 
        // If the scanned character is
        // an operand, add it to output string.
        //     NOTE: This is MOSTLY fine for converting expressions to postfix, but I need to make it 
        // work for numbers bigger than one digit. Possibly require that numeric characters get wrapped
        // in square brackets [like so.] Upon being evaluated, numbers in that format will get pushed 
        // to a stack, and added to each other after multiplying based on the stack depth. For example,
        // '123' would be placed into the output string as [123]. Upon evaluation, it would be read 
        // into a stack as '1', '2', '3', then evaluated numberically as 3*10^0 + 2*10^1 + 1*10^2.
        // ^ BRACKETING PART IS DONE

        if(c >= '0' && c <= '9'){
            if(!s[i-1] || (s[i-1] && !isdigit(s[i-1])) ){ // If the previous character was not a number, open bracket 
                result += '[';
            }

            // Write number
            result += c;

            if (!s[i+1] || !isdigit(s[i+1])){ // If the next character is not a number, close the bracket
                result += ']';
            }
        }
        else if((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))
            result += c;
 
        // If the scanned character is an
        // ‘(‘, push it to the stack.
        else if(c == '(')
            st.push('(');
 
        // If the scanned character is an ‘)’,
        // pop and to output string from the stack
        // until an ‘(‘ is encountered.
        else if(c == ')') {
            while(st.top() != '(')
            {
                result += st.top();
                st.pop();
            }
            st.pop();
        }
 
        //If an equation symbol is scanned
        else if(c == '=' || c == '>' || c == '<'){
            while(!st.empty()) {
                result += st.top();
                st.pop(); 
            }
            result += c;
        }
        else{ // If an operator is scanned
            while(!st.empty() && Prec(s[i]) <= Prec(st.top())) {
                result += st.top();
                st.pop(); 
            }
            st.push(c);
        }
    }
 
    // Pop all the remaining elements from the stack
    while(!st.empty()) {
        result += st.top();
        st.pop();
    }
 
    return result;
}

/**
 * @brief Evaluates a postfix equation for truth
 * @param _eq The equation to evaluate, must be in postfix form and all values must be integers
 * @return True or False depending on if the equation is true or false
 */
bool EvaluatePostEquation(std::string _eq){
    std::stack<int> stack;
    int lhs;
    char symbol;

    // Scan all characters one by one
    for (int i = 0; i < _eq.length(); i++)
    {
        char c = _eq[i];
        //std::cout << "Reading: " << c << "\n";
        // Keep track of where an integer started
        if (c == '['){
            stack.push(c);
        }
        else if (c == ']'){
            // Current sum
            int num = 0;

            // How deep in the stack did I pull the current number from?
            int counter = 0;
            while(stack.top() != '['){
                // Add current number * 10^depth to find actual number value
                num += stack.top() * Pow(10, counter);
                stack.pop();

                // Increment the counter so you know what power to raise the next number to
                counter++;
            }

            // Pop '['
            stack.pop();

            // Push the correct value
            stack.push(num);
        }
        else if (isdigit(c)){
            stack.push(c - '0');
        }
        // If the scanned character is an equation symbol, save
        // left hand side and pop the stack. Save the symbol off
        // the stack.
        else if (c == '=' || c == '<' || c == '>'){
            symbol = c;

            lhs = stack.top();
            stack.pop();
        }
        // If the scanned character is an operator, pop two
        // elements from stack apply the operator
        else
        {   
            // Remember to pop after getting the value
            int val1 = stack.top();
            stack.pop();

            int val2 = stack.top();
            stack.pop();

            switch (_eq[i])
            {
            case '+': stack.push(val2 + val1); break;
            case '-': stack.push(val2 - val1); break;
            case '*': stack.push(val2 * val1); break;
            case '/': stack.push(val2 / val1); break;
            case '^': stack.push(Pow(val2, val1)); break;
            }
        }
    }
    
    //std::cout << "Evaluating equation" << "\n";

    // stack.top() should now be the right hand side of the equation.
    if(symbol == '='){
        return lhs == stack.top();
    }
    else if (symbol == '<'){
        return lhs < stack.top();
    }
    else if (symbol == '>'){
        return lhs > stack.top();
    }

    // return false if the code gets down here somehow
    return false;
}

/**
 * @brief Evaluates a postfix expression for distance from truth, measured in "error"
 * @param _eq The equation to evaluate, must be in postfix form and all values must be integers
 * @return 0 if the equation is true, otherwise abs(lhs - rhs)
 */
int EvaluatePostEquationError(std::string _eq){
    std::stack<int> stack;
    int lhs;
    char symbol;

    // Scan all characters one by one
    for (int i = 0; i < _eq.length(); i++)
    {
        char c = _eq[i];
        // Keep track of where an integer started
        if (c == '['){
            stack.push(c);
        }
        else if (c == ']'){
            // Current sum
            int num = 0;

            // How deep in the stack did I pull the current number from?
            int counter = 0;
            while(stack.top() != '['){
                // Add current number * 10^depth to find actual number value
                num += stack.top() * Pow(10, counter);
                stack.pop();

                // Increment the counter so you know what power to raise the next number to
                counter++;
            }

            // Pop '['
            stack.pop();

            // Push the correct value
            stack.push(num);
        }
        else if (isdigit(c)){
            stack.push(c - '0');
        }
        // If the scanned character is an equation symbol, save
        // left hand side and pop the stack. Save the symbol off
        // the stack.
        else if (c == '=' || c == '<' || c == '>'){
            symbol = c;

            lhs = stack.top();
            stack.pop();
        }
        // If the scanned character is an operator, pop two
        // elements from stack apply the operator
        else
        {   
            // Remember to pop after getting the value
            int val1 = stack.top();
            stack.pop();

            int val2 = stack.top();
            stack.pop();

            // Arithmetic is done with val2 first because pushing and then popping the variables from the stack reversed them
            switch (_eq[i])
            {
            case '+': stack.push(val2 + val1); break;
            case '-': stack.push(val2 - val1); break;
            case '*': stack.push(val2 * val1); break;
            case '/': stack.push(val2 / val1); break;
            case '^': stack.push(Pow(val2, val1)); break;
            }
        }
    }
    

    // stack.top() should now be the right hand side of the equation.
    if(symbol == '='){
        // Return error
        // This will be 0 if lhs and rhs agree
        return abs(lhs - stack.top());
    }
    else if (symbol == '<'){
        if(lhs < stack.top()){
            // No cost if correct
            return 0;
        }
        else{
            // +1 because otherwise it would count as 0 cost if both sides have same value
            return 1 + abs(lhs - stack.top());
        }
    }
    else if (symbol == '>'){
        if(lhs > stack.top()){
            return 0;
        }
        else{
            return 1 + abs(lhs - stack.top());
        }
    }

    // return false if the code gets down here somehow
    return 0;
}

#endif