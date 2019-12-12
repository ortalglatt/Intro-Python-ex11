import random
import itertools


class Node:
    """
    Class of Node objects, contains the data of the node, the negative child
    and the positive child.
    This class contains functions that return information about the node, or
    change the node's attributes.
    """

    def __init__(self, data, pos=None, neg=None):
        """
        Initialize a new Node object.
        :param data: a string that represents the data in the node
        :param pos: None or a Node object - positive child of the node
        :param neg: None or a Node object - negative child of the node
        """
        self.data = data
        self.positive_child = pos
        self.negative_child = neg

    def get_data(self):
        """
        :return: the node's data
        """
        return self.data

    def get_pos_child(self):
        """
        :return: the positive child of the node
        """
        return self.positive_child

    def get_neg_child(self):
        """
        :return: the negative child of the node
        """
        return self.negative_child

    def set_pos_child(self, pos_child):
        """
        :param pos_child: the new positive child
        :return: Change the node's positive child to the given pos_child.
        """
        self.positive_child = pos_child

    def set_neg_child(self, neg_child):
        """
        :param neg_child: the new negative child
        :return: Change the node's negative child to the given neg_child.
        """
        self.negative_child = neg_child


class Record:
    """
    Class of Record objects, illness and a list of symptoms of the record.
    This class contains functions that return information about the record.
    """

    def __init__(self, illness, symptoms):
        """
        Initialize a new Record object.
        :param illness: the record's illness
        :param symptoms: a list of symptoms that lead to the illness
        """
        self.illness = illness
        self.symptoms = symptoms

    def get_illness(self):
        """
        :return: the record's illness
        """
        return self.illness

    def get_symptoms(self):
        """
        :return: the record's symptoms
        """
        return self.symptoms


def parse_data(filepath):
    with open(filepath) as data_file:
        records = []
        for line in data_file:
            words = line.strip().split()
            records.append(Record(words[0], words[1:]))
        return records


class Diagnoser:
    """
    Class of Diagnoser objects, contains the root of the diagnoser tree.
    This class contains functions that return information about the tree.
    """

    HEALTHY = "healthy"

    def __init__(self, root):
        """
        Initialize a new Diagnoser object.
        :param root: a Node object that represent the root of the tree
        """
        self.root = root

    def get_root(self):
        """
        :return: the diagnoser's root
        """
        return self.root

    def diagnose(self, symptoms):
        """
        Diagnoses which illness feats the symptoms according to the tree.
        :param symptoms: a list of symptoms
        :return: the diagnosed illness
        """
        diagnose_illness = self.__diagnose_helper(symptoms, self.root)
        return diagnose_illness

    def __diagnose_helper(self, symptoms, root):
        """
        Recursive function that diagnose the illness that feats the symptoms
        according to the tree.
        :param symptoms: a list of symptoms we want to check.
        :param root: the root of the tree.
        :return: the diagnosed illness.
        """
        if not root.get_pos_child():
            return root.get_data()
        elif root.get_data() in symptoms:
            new_root = root.get_pos_child()
            return self.__diagnose_helper(symptoms, new_root)
        else:
            new_root = root.get_neg_child()
            return self.__diagnose_helper(symptoms, new_root)

    def calculate_success_rate(self, records):
        """
        Calculates the success rate of the tree, by checking if the record's
        symptoms lead the path in the tree to the record's illness.
        :param records: a list of record's objects.
        :return: the success rate of the tree according to the given records.
        """
        success_counter = 0
        for record in records:
            if self.diagnose(record.get_symptoms()) == record.get_illness():
                success_counter += 1
        return success_counter / len(records)

    def all_illnesses(self):
        """
        :return: a sorted list of all the illnesses in the tree (the illness that
        has the biggest amount of appearances in the tree will be first).
        """
        ills_dict = {}
        self.__all_illnesses_helper(ills_dict, self.root)
        illnesses_list = sorted(ills_dict, key=ills_dict.get, reverse=True)
        return illnesses_list

    def __all_illnesses_helper(self, illnesses_dict, root):
        """
        Recursive function that edit the illnesses dictionary, and adds all the
        illnesses in the tree.
        :param illnesses_dict: a dictionary of all the illnesses and their
        amount of appearances in the tree (empty in the beginning).
        :param root: the root af the tree.
        :return: None
        """
        if not root.get_pos_child():
            if root.get_data() in illnesses_dict:
                illnesses_dict[root.get_data()] += 1
            elif root.get_data() != self.HEALTHY:
                illnesses_dict[root.get_data()] = 1
            return
        new_root = root.get_pos_child()
        self.__all_illnesses_helper(illnesses_dict, new_root)
        new_root = root.get_neg_child()
        self.__all_illnesses_helper(illnesses_dict, new_root)

    def most_rare_illness(self, records):
        """
        Checks which illness in the tree is the most rare illness, by checking
        how many records lead to each illness.
        :param records: a list of record's object.
        :return: the most rare illness.
        """
        ills_dic = {illness: 0 for illness in self.all_illnesses()}
        for record in records:
            diagnose_illness = self.diagnose(record.get_symptoms())
            if diagnose_illness != self.HEALTHY:
                ills_dic[diagnose_illness] += 1
        return min(ills_dic, key=ills_dic.get)

    def paths_to_illness(self, illness):
        """
        :param illness: the illness we want to get all the paths that gets to
        it in the tree.
        :return: a list of all the paths that will bring you to the illness.
        """
        paths_lst = []
        path = []
        self.__paths_to_illness_helper(self.root, illness, paths_lst, path)
        return paths_lst

    def __paths_to_illness_helper(self, root, illness, paths_lst, path):
        """
        Recursive function that adds all the paths to the given illness, to the
        given paths list.
        :param root: the root of the tree.
        :param illness: the illness we want to find all the paths to.
        :param paths_lst: a list of all the paths to the illness (empty on the
        beginning).
        :param path: a list of "True"s and "False" that symbolize the current
        path.
        :return: None
        """
        if not root.get_pos_child():
            if root.get_data() == illness:
                paths_lst.append(path)
                return
            return
        new_root = root.get_pos_child()
        self.__paths_to_illness_helper(new_root, illness, paths_lst,
                                                      path + [True])
        new_root = root.get_neg_child()
        self.__paths_to_illness_helper(new_root, illness, paths_lst,
                                                      path + [False])


def build_tree(records, symptoms):
    """
    Builds a tree with all the symptoms, that will bring you to the illness
    that you have the most chances to have.
    :param records: a list of record's objects
    :param symptoms: a list of symptoms
    :return: the root of the tree we built.
    """
    root = Node(symptoms[0])
    path_to_node = []
    build_tree_helper(records, symptoms, root, path_to_node)
    return root


def build_tree_helper(records, symptoms, node, path_to_node, index=1):
    """
    Recursive function that builds the tree with all the given symptoms, that
    leads the user to the illness according to his path and to the given
    records.
    :param records: a alis of record's objects
    :param symptoms: a list of symptoms
    :param node: the node the tree wil start from
    :param path_to_node: the current path until the node
    :param index: at the beginning 1
    :return: None
    """
    if index == len(symptoms):
        pos_data = choose_illness(records, path_to_node +
                                  [(True, symptoms[index - 1])])
        node.set_pos_child(Node(pos_data))
        neg_data = choose_illness(records, path_to_node +
                                  [(False, symptoms[index - 1])])
        node.set_neg_child(Node(neg_data))
        return
    pos_node = Node(symptoms[index])
    neg_node = Node(symptoms[index])
    node.set_pos_child(pos_node)
    node.set_neg_child(neg_node)
    build_tree_helper(records, symptoms, pos_node, path_to_node +
                      [(True, symptoms[index - 1])], index + 1)
    build_tree_helper(records, symptoms, neg_node, path_to_node +
                      [(False, symptoms[index - 1])], index + 1)


def choose_illness(records, path):
    """
    Find the illness the feats most to the path according to the given records.
    If their is no illness that feats, the function will choose a random
    illness from the record's illnesses.
    :param records: a list of record's objects
    :param path: a list of tuples of all the symptoms with True or False,
    according to the way on the tree.
    :return: illness that feats the path and the records.
    """
    ills_dic = {}
    for record in records:
        if does_record_feat(record, path):
            if record.get_illness() in ills_dic:
                ills_dic[record.get_illness()] += 1
            else:
                ills_dic[record.get_illness()] = 1
    if ills_dic:
        return max(ills_dic, key=ills_dic.get)
    all_illnesses = [record.get_illness() for record in records]
    return random.choice(all_illnesses)


def does_record_feat(record, path):
    """
    Check if the symptoms according to the path feats the symptoms in the given
    record.
    :param record: the record we want to check
    :param path: a list of tuple of all the symptoms with True or False,
    according to the way on the tree.
    :return: True if the path feats the record's symptoms, and False if it's
    not.
    """
    for step in path:
        if step[0]:
            if step[1] not in record.get_symptoms():
                return False
        else:
            if step[1] in record.get_symptoms():
                return False
    return True


def optimal_tree(records, symptoms, depth):
    """
    Builds trees of all the symptoms' sets that their size is the given depth.
    The function checks which tree has the highest success rate.
    :param records: a list of record's objects
    :param symptoms: a list of symptoms
    :param depth: int, the depth of the tree we want to check
    :return: the root of the optimal tree
    """
    all_comb = [comb for comb in itertools.combinations(symptoms, depth)]
    success_dic = {}
    for comb in all_comb:
        diag = Diagnoser(build_tree(records, comb))
        success_dic[diag] = diag.calculate_success_rate(records)
    return max(success_dic, key=success_dic.get).get_root()

