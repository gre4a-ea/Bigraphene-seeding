#include <iostream>


class Product{
private:
  std::string name;
  double cost;
  double weight;
public:
  Product(std::string name, double cost, double weight){
    this -> name = name;
    this -> cost = cost;
    this -> weight = weight;
  }

  virtual std::string getName() { return name; }

  virtual double getCost() { return cost; }

  virtual double getWeight() { return weight; }

  virtual void setName(std::string n) { this -> name = n; }

  virtual void setCost_pr(double c) { this -> cost = c; }

  virtual void setWeight_pr(double w) { this -> weight = w; }

};


class Buy: public virtual Product{
private:
  int number_of_products;
  double total_cost;
  double total_weight;
public:
  Buy(std::string name,
    double cost,
    double weight,
    int number_of_products): Product(name, cost, weight){
    this -> number_of_products = number_of_products;
    this -> total_cost = number_of_products * cost;
    this -> total_weight = number_of_products * weight;
  }

  virtual void setNumber_of_products(int n) {
    this -> number_of_products = n;
    setGetTotal_cost();
    setGetTotal_weight();
  }

  virtual void setGetTotal_cost() {
    this -> total_cost = number_of_products * Product::getCost();
  }

  virtual void setGetTotal_weight() {
    this -> total_weight = number_of_products * Product::getWeight();
  }

  virtual int get_Number_of_products() { return number_of_products; }

  virtual double get_Total_cost() { return total_cost; }

  virtual double get_Total_weight() { return total_weight; }
};


class Check: public virtual Buy, public virtual Product{
public:
  Check(std::string name,
    double cost,
    double weight,
    int number_of_products):
  Buy(name, cost, weight, number_of_products),
  Product(name, cost, weight)
  {}

  void print_Name() {
    std::cout << "Good: " << Product::getName() << '\n';
  }

  void print_Cost() {
    std::cout << "Cost: " << Product::getCost() << " rub" << '\n';
  }

  void print_Weight() {
    std::cout << "Weight: " << Product::getWeight() << " g" << '\n';
  }

  void print_Number_of_products() {
    std::cout << "Quantity of goods: " << Buy::getNumber_of_products() << " grand" << '\n';
  }

  void print_Total_cost() {
    std::cout << "Total Cost: " << Buy::getTotal_cost() << " rub" << '\n';
  }

  void print_Total_weight() {
    std::cout << "Total Weight: " << Buy::getTotal_weight() << " g" << '\n';
  }

  void print_AllAboutBuy() {
    print_Number_of_products();
    print_Total_cost();
    print_Total_weight();
  }

  void print_AllAboutProduct() {
    print_Name();
    print_Cost();
    print_Weight();
  }

  void setCost(double c) {
    Product::setCost_pr(c);
    Buy::setGetTotal_cost();
  }

  void setWeight(double w) {
    Product::setWeight_pr(w);
    Buy::setGetTotal_weight();
  }
};


int main() {
  std::string s = "Book";

  Check checker(s, 1000, 600, 5);

  checker.print_AllAboutProduct();
  std::cout << '\n';
  checker.print_AllAboutBuy();

  std::cout << '\n';

  checker.setNumber_of_products(10);
  checker.print_AllAboutProduct();
  std::cout << '\n';
  checker.print_AllAboutBuy();

  std::cout << '\n';

  checker.setWeight(200);
  checker.print_AllAboutProduct();
  std::cout << '\n';
  checker.print_AllAboutBuy();
  return 0;
}
