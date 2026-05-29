#imports
import math
import pytest
from Homework_6_vector_class import Vector

#tests
def test_eq_and_ineq():
    v = Vector(1, 2, 3)
    w = Vector(1, 2, 3)
    u = Vector(3, 2, 1)

    assert v == w
    assert v != u


def test_add():
    v = Vector(1, 2, 3)
    w = Vector(4, 5, 6)
    assert v + w == Vector(5, 7, 9)


def test_sub():
    v = Vector(5, 5, 5)
    w = Vector(1, 2, 3)
    assert v - w == Vector(4, 3, 2)


def test_dot():
    v = Vector(1, 2, 3)
    w = Vector(4, 5, 6)
    assert v * w == 32


def test_cross():
    v = Vector(1, 0, 0)
    w = Vector(0, 1, 0)
    assert v.cross(w) == Vector(0, 0, 1)


def test_length():
    v = Vector(3, 4, 12)
    assert v.length() == math.sqrt(169)


def test_hash_and_set():
    v = Vector(1, 2, 3)
    w = Vector(1, 2, 3)
    u = Vector(4, 5, 6)

    s = {v, w, u}
    assert len(s) == 2