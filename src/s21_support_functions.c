#include "s21_matrix.h"

int s21_is_not_correct_matrix(matrix_t *result) {
  ResultVal status = EXIT_OK;
  if (!result || !result->matrix || result->columns < 1 || result->rows < 1)
    status = EXIT_INCORR_MATRIX;

  return status;
}

int s21_has_equal_sizes(matrix_t *A, matrix_t *B) {
  return ((A->columns == B->columns) && (A->rows == B->rows) ? 1 : 0);
}

double s21_calc_det(matrix_t *A) {
  double result = 0;
  if (A->rows == 1) {
    result = A->matrix[0][0];
  } else if (A->rows == 2) {
    result = (A->matrix[0][0] * A->matrix[1][1]) -
             (A->matrix[0][1] * A->matrix[1][0]);
  } else {
    for (int j = 0; j < A->columns; ++j) {
      matrix_t tmp;
      s21_create_matrix(A->rows - 1, A->columns - 1, &tmp);
      s21_get_minor(A, 0, j, &tmp);
      if (!(j & 1)) {
        result += A->matrix[0][j] * s21_calc_det(&tmp);
      } else {
        result -= A->matrix[0][j] * s21_calc_det(&tmp);
      }
      s21_remove_matrix(&tmp);
    }
  }

  return result;
}

// Функция для вычисления минора с учетов номера столбца и строки
void s21_get_minor(matrix_t *A, int i, int j, matrix_t *tmp) {
  int tmp_row = 0;
  int tmp_col = 0;
  for (int row = 0; row < A->rows; ++row) {
    for (int col = 0; col < A->columns; ++col) {
      if (col != j && row != i) {
        tmp->matrix[tmp_row][tmp_col] = A->matrix[row][col];
        ++tmp_col;
      }
    }
    if (row != i) {
      ++tmp_row;
      tmp_col = 0;
    }
  }
}

// Вычисление матрицы миноров
void s21_calc_minor(matrix_t *A, matrix_t *result) {
  for (int i = 0; i < A->rows; ++i) {
    for (int j = 0; j < A->columns; ++j) {
      matrix_t tmp;
      s21_create_matrix(A->rows - 1, A->columns - 1, &tmp);
      s21_get_minor(A, i, j, &tmp);
      result->matrix[i][j] = s21_calc_det(&tmp);
      s21_remove_matrix(&tmp);
    }
  }
}