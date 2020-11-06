# LL(1) parser implementation for the module definition of KEGGs.
#
# It ultimately allows, after having parsed a definition, to check 
# how many kegg orthologs are missing to have a complete module.
#
# For now, it does not support signature modules.
#
# TODO (BM): classes naming is confusing, should have clear difference between
# tokens and nodes of the syntax tree.
# 
# Author: bastian.marquis@protonmail.com
# Date: 05.11.2020

class KoModuleTree:
    def __init__(self, root):
        self.root = root

    def tokenize():
        pass

    def parse(ko_definition):

        module_def = KoModule(root)
        return module_def

    def get_n_missing(self, kos):
        return self.root.get_n_missing(kos)

    def get_ko_ids(self):
        return self.get_ko_ids()

class Complex:
    def __init__(self, list_comp):
        self.list_comp = list_comp

    def get_n_missing(self, kos):
        return sum(node.get_n_missing() for node in self.list_comp)

    def get_ko_ids(self):
        return (node.get_ko_ids() for node in list_comp)

class Component:
    def __init__(self, ko_id):
        self.ko_id = ko_id

    def get_n_missing(self, kos):
        if ko_id not in kos:
            return 1
        else:
            return 0

    def get_ko_ids(self):
        return self.ko_id

class OptionalSubunit(Complex):
    def get_n_missing(self, kos):
        return 0

# May need to put those in a different file
class KoAnd:
    def __init__(self, list_and):
        self.list_and = list_and

    def get_n_missing(self, kos):
        n_missing = 0
        for node in self.list_and:
            n_missing += node.get_n_missing(kos)
        return n_missing

    def get_ko_ids(self):
        return (node.get_ko_ids() for node in self.list_and)

class KoOr:
    def __init__(self, list_or):
        self.list_or = list_or

    def get_n_missing(kos):
        return min(node.get_n_missing() for node in self.list_or)

    def get_ko_ids(self):
        return (node.get_ko_ids() for node in self.list_and)

class KoNode:
    def __init__(self, node_id):
        self.node_id = node_id

    def get_n_missing(self, kos):
        if self.node_id in kos:
            return 0
        return 1

    def get_ko_ids(self):
        return self.node_id

class UndefinedKoNode:
    def get_n_missing(self):
        return 0

    def get_ko_ids(self):
        return []

class Token:
    pass

class LeftParToken(Token):
    def __str__(self):
        return "LeftPar"

class RightParToken(Token):
    def __str__(self):
        return "RightPar"

class Ko(Token):
    def __init__(self, ko_id):
        self.ko_id = ko_id

    def __str__(self):
        return f"KO{self.ko_id}"

class AndToken(Token):
    def __str__(self):
        return "And"

class ComplexToken(Token):
    def __str__(self):
        return "Complex"

class OptionalComplexComponent(Token):
    def __str__(self):
        return "OptComplex"

class OrToken(Token):
    def __str__(self):
        return "Or"

class UndefinedToken(Token):
    def __str__(self):
        return "Undefined KO"


class Tokenizer:
    KO_TOKEN_LENGTH = 5

    def __init__(self, module_def):
        self.module_def = module_def

    def __iter__(self):
        self.i = 0
        self.ttl_len = len(self.module_def)
        return self

    def tokenize_ko(self):
        self.i += 1
        acc = 0
        for i in range(Tokenizer.KO_TOKEN_LENGTH):
            curr_char = self.module_def[self.i]
            if not curr_char.isdigit():
                raise Exception("Invalid KO")
            acc = acc*10 + int(curr_char)
            self.i += 1
        return Ko(acc)

    def tokenize_minus(self):
        self.i += 1
        if self.module_def[self.i] == "-":
            self.i += 1
            return UndefinedToken()
        elif self.module_def[self.i] == "K":
            return OptionalComplexComponent()
        elif self.module_def[self.i] == "(":
            return OptionalComplexComponent()
        else:
            raise Exception("Error tokenizing a minus")

    def __next__(self):
        while self.i < self.ttl_len:
            char = self.module_def[self.i]
            if char == "(":
                self.i += 1
                return LeftParToken()
            elif char == ")":
                self.i += 1
                return RightParToken()
            elif char == " ":
                self.i += 1
                return AndToken()
            elif char == "+":
                self.i += 1
                return ComplexToken()
            elif char == ",":
                self.i += 1
                return OrToken()
            elif char == "-":
                return self.tokenize_minus()
            elif char == "K":
                return self.tokenize_ko()
            else:
                raise Exception(f"Error tokenizing at position {self.i}({self.module_def[self.i]})")
        raise StopIteration


class ModuleParser:

    def __init__(self, module_def):
        self.module_def = module_def
        tokenizer = Tokenizer(module_def)
        self.token_iter = iter(tokenizer)
        self.curr_token = None

    def next_token(self):
        if self.curr_token != None:
            token, self.curr_token = self.curr_token, None
            return token
        try:
            return next(self.token_iter)
        except StopIteration as e:
            return None

    def parse_parentheses(self):
        subtree = self.parse()
        next_token = self.next_token()
        if next_token == None or not isinstance(next_token, RightParToken):
            raise Exception("Unexpected end of expression (unbalanced parentheses) " + str(next_token))
        return subtree

    def parse_leaf(self):
        n = self.next_token()

        if isinstance(n, Ko):
            return KoNode(n.ko_id)
        elif isinstance(n, UndefinedToken):
            return UndefinedKoNode()
        elif isinstance(n, LeftParToken):
            return self.parse_parentheses()
        elif isinstance(n, OptionalComplexComponent):
            return OptionalSubunit([self.parse_leaf()])
        else:
            raise Exception("Unexpected token: "+str(n))

    def parse_complex(self):
        node = self.parse_leaf()
        n = self.next_token()

        if isinstance(n, ComplexToken):
            compl = Complex([node])
            compl.list_comp.append(self.parse())
            return compl
        elif isinstance(n, OptionalComplexComponent):
            optional_compl = OptionalSubunit([node])
            optional_compl.list_comp.append(self.parse())
            return optional_compl
        else:
            self.curr_token = n
            return node

    def parse(self):
        comp = self.parse_complex()
        n = self.next_token()
        if n==None:
            return comp
        elif isinstance(n, AndToken):
            and_node = KoAnd([comp])
            and_node.list_and.append(self.parse())
            return and_node
        elif isinstance(n, OrToken):
            or_node = KoOr([comp])
            or_node.list_or.append(self.parse())
            return or_node
        else:
            self.curr_token = n
            return comp
