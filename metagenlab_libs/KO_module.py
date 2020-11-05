

class KoModuleTree:
    def KoModuleTree(self, root):
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


# May need to put those in a different file
class KoAnd:
    def KoAnd(self, list_and):
        self.list_and = list_and

    def get_n_missing(self, kos):
        n_missing = 0
        for node in self.list_and:
            n_missing += node.get_n_missing(kos)
        return n_missing

    def get_ko_ids(self):
        return (node.get_ko_ids() for node in self.list_and)

class KoOr:
    def KoOr(self, list_or):
        self.list_or = list_or

    def get_n_missing(kos):
        return min(node.get_n_missing() for node in self.list_or)

    def get_ko_ids(self):
        return (node.get_ko_ids() for node in self.list_and)

class KoNode:
    def KoNode(self, node_id):
        self.node_id = node_id

    def get_n_missing(self, kos):
        if self.node_id in kos:
            return 0
        return 1

    def get_ko_ids(self):
        return self.node_id

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