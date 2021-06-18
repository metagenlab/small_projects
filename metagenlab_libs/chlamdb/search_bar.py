
import os

from collections import namedtuple

from whoosh.fields import SchemaClass, TEXT, KEYWORD, ID
from whoosh.qparser import MultifieldParser
from whoosh import index


class SearchBarSchema(SchemaClass):
    locus_tag = KEYWORD(stored=True)
    gene = KEYWORD(stored=True)
    product = TEXT(stored=True)
    organism = TEXT(stored=True)
    cog = KEYWORD(stored=True)
    ko = KEYWORD(stored=True)
    og = KEYWORD(stored=True)
    pfam = KEYWORD(stored=True)


SearchResult = namedtuple("SearchResult", ChlamdbIndex.Field_list)


class ChlamdbIndex:
    Field_list = ["locus_tag", "gene", "product", "organism", "cog", "ko", "og", "pfam"]

    def new_index(name):
        chlamdb_index = ChlamdbIndex()
        chlamdb_index.index = index.create_in(name, SearchBarSchema)
        chlamdb_index.writer = chlamdb_index.index.writer()
        return chlamdb_index

    def use_index(name):
        chlamdb_index = ChlamdbIndex()
        chlamdb_index.index = index.open_dir(name)
        return chlamdb_index

    def add(self, **kwargs):
        self.writer.add_document(**kwargs)

    def search(self, user_query, limit=10):
        parser = MultifieldParser(ChlamdbIndex.Field_list, self.index.schema)
        query = parser.parse(user_query)

        for result in self.index.searcher().search(query, limit=limit):
            gene = result.get("gene", None)
            locus_tag = result.get("locus_tag", None)
            product = result.get("product", None)
            organism = result.get("organism", None)
            cog = result.get("cog", None)
            ko = result.get("ko", None)
            og = result.get("og", None)
            pfam = result.get("pfam", None)
            yield SearchResult(locus_tag=locus_tag, gene=gene, product=product,
                    organism=organism, cog=cog, ko=ko, og=og, pfam=pfam)

    def done_adding(self):
        self.writer.commit(optimize=True)
